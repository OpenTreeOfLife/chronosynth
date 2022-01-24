"""Functions for working with chronograms from Phylesystem"""
#!/usr/bin/env python3
import sys
import os
import opentree
import datetime
import json
import subprocess
import re
import dendropy
import random
import configparser
from sh import git

from opentree import OT

import logging

import chronosynth

from peyotl.phylesystem.git_actions import PhylesystemGitAction

import peyotl

config = configparser.ConfigParser()
config.read(chronosynth.configfile)


log = logging.getLogger(__name__)
log.debug("logging set to debug")


DC = opentree.object_conversion.DendropyConvert()


def set_dev():
    global OT
    OT = opentree.ot_object.OpenTree(api_endpoint='dev')


def set_prod():
    global OT
    OT = opentree.ot_object.OpenTree()

def print_endpoint():
    print(OT._api_endpoint)
    log.debug("api_endpoint is %s", format(OT._api_endpoint))

def find_trees(search_property="ot:branchLengthMode", value="ot:time"):
    """
    Get study and tree ids for all chronograms (trees with branch lengths proportional
    to time) in Phylesystem, i.e., where 'ot:branchLengthMode' == 'ot:time'

    Returns
    -------
    A list of Phylesystem chronogram's source ids: study_id@treeid.
    """
    output = OT.find_trees(search_property=search_property, value=value)
    chronograms = set()
    for study in output.response_dict["matched_studies"]:
        study_id = study['ot:studyId']
        for tree in study['matched_trees']:
            tree_id = tree['ot:treeId']
            chronograms.add('{}@{}'.format(study_id, tree_id))
    return list(chronograms)

def as_dendropy(source_id):
    """
    Get a dendropy object of a chronogram in Phylesystem.

    Example
    -------
    source_id = 'ot_1000@tree1'
    """
    assert '@' in source_id
    study_id = source_id.split('@')[0]
    tree_id = source_id.split('@')[1]
    study = OT.get_study(study_id)## Todo: catch failure of study GET
    study_nexson = study.response_dict['data']
    tree_obj = DC.tree_from_nexson(study_nexson, tree_id)
    return tree_obj

def node_ages(source_id, ultrametricity_precision=None):
    """
    Get node ages for a DendroPy tree object.

    Returns
    -------
    A dictionary of dictionaries with
    {'metadata': {'study_id': study_id, 'tree_id': tree_id, 'time_unit': time_unit},
    'ages': {node_label:node_age}
    }

    """
    if ultrametricity_precision == None:
        ultrametricity_precision = float(config.get('params', 'ultrametricity_precision',
                                     fallback='0.01'))

    study_id = source_id.split('@')[0]
    tree_id = source_id.split('@')[1]
    dp_tree = as_dendropy(source_id)
    assert dp_tree.annotations.get_value("branchLengthMode") == 'ot:time'
    time_unit = dp_tree.annotations.get_value("branchLengthTimeUnit")
    metadata = {'study_id': study_id, 'tree_id': tree_id, 'time_unit' :time_unit}
    ages = {}
    dp_tree.internal_node_ages(ultrametricity_precision=ultrametricity_precision)
    for node in dp_tree.internal_nodes():
        assert node.label not in ages.keys()
        ages[node.label] = node.age
    ret = {'metadata':metadata, 'ages':ages}
    return ret


def map_conflict_ages(source_id, ultrametricity_precision=None):
    """
    Takes a source id in format study_id@tree_id

    returns a dictionary of:
    {'metadata':{'study_id': study_id, 'tree_id': tree_id,
                 'time_unit': time_unit, 'synth_version':version},
    'supported_nodes':{synth_node_id : {'age':age, 'node_label':node_label}}
    }
    """
    ages_data = node_ages(source_id, ultrametricity_precision=ultrametricity_precision)
    metadata = ages_data['metadata']
    version = OT.about()
    metadata['synth_tree_about'] = version['synth_tree_about']
    output_conflict = OT.conflict_info(study_id=metadata['study_id'], tree_id=metadata['tree_id'])
    conf = output_conflict.response_dict
    if(conf==None):
        url = "https://tree.opentreeoflife.org/curator/study/view/{}/?tab=home&tree={}".format(metadata['study_id'], metadata['tree_id'])
        sys.stderr.write("No conflict data available for tree {} \n Check its status at {}\n".format(source_id, url))
        return(None)
    supported_nodes = {}
    for node_label in ages_data['ages']:
        age = ages_data['ages'][node_label]
        if node_label not in conf:
            # This not only happens for the root
            # TODO: map root to synth using mrca??
            ## Skips not in ingroup...
            # print(node_label)
            continue
        node_conf = conf[node_label]
        status = node_conf['status']
        witness = node_conf['witness']
        if status == 'supported_by':
            supported_nodes[witness] = {'age':age, 'node_label':node_label}
    ret = {'metadata':metadata, 'supported_nodes':supported_nodes}
    return ret


def map_conflict_nodes(source_id):
    """
    Takes a source id in format study_id@tree_id

    returns a dictionary of:
    {
    'matched_nodes':{node_id : synnode_label}}
    }
    """
    study_id = source_id.split('@')[0]
    tree_id = source_id.split('@')[1]
    dp_tree = as_dendropy(source_id)
    time_unit = dp_tree.annotations.get_value("branchLengthTimeUnit")
    metadata = {'study_id': study_id, 'tree_id': tree_id}
    output_conflict = OT.conflict_info(study_id, tree_id)
    conf = output_conflict.response_dict
    if(conf==None):
        url = "https://tree.opentreeoflife.org/curator/study/view/{}/?tab=home&tree={}".format(metadata['study_id'], metadata['tree_id'])
        sys.stderr.write("No conflict data available for tree {} \n Check its status at {}\n".format(source_id, url))
        return(None)
    matched_nodes = {}
    for node in dp_tree:
        if node.label in conf:
            node_conf = conf[node.label]
            status = node_conf['status']
            witness = node_conf['witness']
            if status == 'supported_by':
                matched_nodes[node.label] = witness
    ret = {'metadata':metadata, 'matched_nodes':matched_nodes, 'tree':dp_tree}
    return ret



def combine_ages_from_sources(source_ids, ultrametricity_precision=None, json_out = None):
#temporary name
    
    """
    inputs: a list of source ids in format study_id@tree_id

    Outputs: a dictionary
    Key: node_id (or ott id) in synth tree
    example
     synth_node_ages['metadata']{'synth_info':12.3, 
                            'date_processed':Today}
     synth_node_ages['node_ages']:
     {mrcaott123ott456 : {
                          [
                          {'source_id':"ot_1000@tree1", 
                          'age': 48, 
                          'node_label':node_label
                           ],
                           }

    """

    synth_node_ages = {'metadata':{}, 'node_ages':{}}

    versions = OT.about()
    synth_node_ages['metadata']['synth_tree_about'] = versions['synth_tree_about']
    synth_node_ages['metadata']['date'] = str(datetime.date.today())
    synth_node_ages['metadata']['phylesystem_sha'] = get_phylesystem_sha()
    for tag in source_ids:
        try:
            res = map_conflict_ages(tag, ultrametricity_precision=ultrametricity_precision)
        except ValueError:
            time_unit = res['metadata']['time_unit']
            log.info('{}, conflict error, {}\n'.format(tag, time_unit))
            continue
        if res==None:
            log.info('{}, conflict empty\n'.format(tag))
        else:
            sys.stdout.write("study {} has {} supported nodes\n".format(tag, len(res["supported_nodes"])))
            source_id = tag
            # assert synth_node_ages['metadata']['synth_tree_about'] == res['metadata']['synth_tree_about']
            time_unit = res['metadata']['time_unit']
            if time_unit == 'Myr':
                assert tag == "{}@{}".format(res['metadata']['study_id'],res['metadata']['tree_id'])
                for synth_node in res['supported_nodes']:
                    if synth_node not in synth_node_ages['node_ages']: #if not a record yet, we need to create one
                        synth_node_ages['node_ages'][synth_node] = []
                    age = res['supported_nodes'][synth_node]['age']
                    source_node = res['supported_nodes'][synth_node]['node_label']
                    entry = {'source_id': source_id, 'age':age, 'source_node':source_node}
                    synth_node_ages['node_ages'][synth_node].append(entry)
            else:
                #skips all tree not in mya
                pass

    sf = json.dumps(synth_node_ages, sort_keys=True, indent=2, separators=(',', ': '), ensure_ascii=True)
    if json_out is not None:
        ofi = open(json_out,'w')
        ofi.write(sf)
        ofi.close()
    return(synth_node_ages)


def get_trees_from_nexson():

def update_ages_from_sources(dates, source_ids, new_sha, ultrametricity_precision=None, json_out = None):
    dates['synth_node_ages']['metadata']['phylesystem_sha'] = new_sha

    for tag in source_ids:
        for synth_node in dates['synth_node_ages']:
            if tag in  dates['synth_node_ages'][synth_node]:
                print(dates['synth_node_ages'][synth_node])
        try:
            res = map_conflict_ages(tag, ultrametricity_precision=ultrametricity_precision)
        except ValueError:
            time_unit = res['metadata']['time_unit']
            log.info('{}, conflict error, {}\n'.format(tag, time_unit))
            continue
        if res==None:
            log.info('{}, conflict empty\n'.format(tag))
        else:
            sys.stdout.write("study {} has {} supported nodes\n".format(tag, len(res["supported_nodes"])))
            source_id = tag
            # assert synth_node_ages['metadata']['synth_tree_about'] == res['metadata']['synth_tree_about']
            time_unit = res['metadata']['time_unit']
            if time_unit == 'Myr':
                assert tag == "{}@{}".format(res['metadata']['study_id'],res['metadata']['tree_id'])
                for synth_node in res['supported_nodes']:
                    if synth_node not in synth_node_ages['node_ages']: #if not a record yet, we need to create one
                        synth_node_ages['node_ages'][synth_node] = []
                    age = res['supported_nodes'][synth_node]['age']
                    source_node = res['supported_nodes'][synth_node]['node_label']
                    entry = {'source_id': source_id, 'age':age, 'source_node':source_node}
                    synth_node_ages['node_ages'][synth_node].append(entry)
            else:
                #skips all tree not in mya
                pass
    sf = json.dumps(synth_node_ages, sort_keys=True, indent=2, separators=(',', ': '), ensure_ascii=True)
    if json_out is not None:
        ofi = open(json_out,'w')
        ofi.write(sf)
        ofi.close()
    return(synth_node_ages)




#This should probably go in peyotl or somethings
def get_phylesystem_sha(repo_url = "https://github.com/OpenTreeOfLife/phylesystem-1.git", repo_dir=None):
    if repo_dir:
        assert(os.path.exists(repo_dir))
        sha = peyotl.git_storage.git_action.get_HEAD_SHA1('{}/.git'.format(repo_dir))
    else:
        process = subprocess.Popen(["git", "ls-remote", repo_url], stdout=subprocess.PIPE)
        stdout, stderr = process.communicate()
        sha = re.split(r'\t+', stdout.decode('ascii'))[0]
    return sha


def pull_phylesystem(repo_dir, repo_url = "https://github.com/OpenTreeOfLife/phylesystem-1.git"):
    assert(os.path.exists(repo_dir))
    git_dir_arg = "--git-dir={}/.git".format(repo_dir)
    git(git_dir_arg, 'pull', repo_url)


def build_synth_node_source_ages(cache_file_path=None, remote=True, ultrametricity_precision=None, repo_dir=None):
    if cache_file_path == None:
        cache_file_path = config.get('paths', 'cache_file_path',
                                     fallback='/tmp/node_ages.json')
    if os.path.exists(cache_file_path):
        dates = json.load(open(cache_file_path))
        current_sha = get_phylesystem_sha(repo_dir=repo_dir)
        # always remote?? 
        # find trees relies on otindex, which relies on github... so maybe remote is best even if a lil slo?
        cached_sha = dates['metadata'].get('phylesystem_sha')
        if cached_sha == current_sha:
            sys.stdout.write("No new changes to phylesystem, using cached dates at {}\n".format(cache_file_path))
            return dates
        if cached_sha != current_sha:
            if repo_dir:
                repo = PhylesystemGitAction(repo=repo_dir)
                changed_studies = repo.get_changed_docs(cached_sha)


            else:
                sys.stdout.write("Phylesystem has changes since dates were cached, reloading and saving to {}\n".format(cache_file_path))
                sources = find_trees()
                dates = combine_ages_from_sources(sources, ultrametricity_precision=ultrametricity_precision, json_out = cache_file_path)
                return dates
    else:
        sources = find_trees()
        dates = combine_ages_from_sources(sources, json_out = cache_file_path)
    return dates

def synth_node_source_ages(node, cache_file_path=None):
    if cache_file_path == None:
        cache_file_path = config.get('Paths', 'cache_file_path',
                                     fallback='/tmp/node_ages.json')
    log.debug("cache file path %s" % cache_file_path)
    ##check if node is in synth?
    synth_resp = OT.synth_node_info(node)
    retdict = {}
    retdict['query'] = node
    if synth_resp.status_code == 200:
        resp_node = synth_resp.response_dict[0]['node_id']
        retdict['synth_node_id'] = resp_node
        if node.startswith('ott'):
            if node != resp_node:
                msg = "Taxon {} is not monophyletic.\
                    resolving to MRCA: synth_node {}, and reporting dates for that\n".format(node, synth_resp)
                retdict['msg'] =  msg
        dates = build_synth_node_source_ages(cache_file_path)
        retdict['ot:source_node_ages'] = dates['node_ages'].get(node)
    elif node.startswith('ott') and node.strip('ott').isnumeric():
        tax_resp = OT.taxon_info(node)
        if tax_resp.status_code == 200:
            msg = "Taxon {} is in the taxonomy, but cannot be found in the synth tree\n".format(node)
            retdict['msg'] =  msg
            retdict = {'msg': msg, 'synth_response': synth_resp.response_dict, 'tax_response': tax_resp.response_dict }
        else:
            msg = "node {} not found in synthetic tree or taxonomy".format(node)
            retdict = {'msg': msg, 'synth_response': synth_resp.response_dict, 'tax_response': tax_resp.response_dict}
    else:
            msg = "node {} not found in synthetic tree or taxonomy.".format(node)
            retdict = {'msg': msg, 'synth_response': synth_resp.response_dict}
    return retdict


def write_fastdate_prior(subtree, dates, var_mult = 0.1, outputfile='node_prior.txt'):
    mrca_dict = {}
    for node in subtree:
        lab = None
        if node.label:
            lab = str(node.label)
            if lab in dates['node_ages']:
                mrca_dict[lab] = {}
                mrca_dict[lab]['ages'] = dates['node_ages'][lab]
                nd = subtree.find_node_with_label(lab)
                mrca_dict[lab]['tips'] = [ti.taxon.label.replace(' ','_') for ti in nd.leaf_iter()]
    if len(mrca_dict) == 0:
        sys.stderr.write("no calibrations")
        return None
    fi=open(outputfile,'w')
    for dated_node in mrca_dict:
        fi.write("'")
        fi.write("','".join(mrca_dict[dated_node]['tips']))
        fi.write("'")
        fi.write(' ')
        ages = [source['age'] for source in mrca_dict[dated_node]['ages']]
        avgage = sum(ages)/len(ages)
        if len(ages) > 1:
            var = statistics.variance(ages)
        else:
            var = var_mult*avgage
        fi.write('norm({},{},{})\n'.format(0, var, avgage))
    fi.close()
    return outputfile

def write_fastdate_tree(subtreepath, br_len = 0.01, polytomy_br = 0.001, outputfile='fastdate_input.tre'):
    subtree=dendropy.Tree.get_from_path(subtreepath, schema="newick")
    subtree.resolve_polytomies(rng=random)    
    subtree.suppress_unifurcations()
    for edge in subtree.levelorder_edge_iter():
        if (edge.tail_node is not None) and (edge.length is None):
            edge.length = 0.01
        if edge.length == 0:
            edge.length = 0.001
    subtree.write(path=outputfile, schema="newick")

def date_synth_subtree(ott_id, reps, max_age=None, summary='sumtre.tre'):
    dates = build_synth_node_source_ages(ultrametricity_precision=0.001)

    if max_age:
            max_age = max_age
    elif 'ott'+ott_id in dates['node_ages']:
        max_age = max([source['age'] for source in dates['node_ages']['ott'+ott_id]]) * 1.25
    else:
        sys.stderr.write("ERROR: no age estimate for root - please provide max root age using --max_age")
        return None

    output = OT.synth_subtree(ott_id=ott_id, label_format='id')
    subtree = dendropy.Tree.get_from_string(output.response_dict['newick'], schema = 'newick')
    subtree.write(path="unresolved.tre", schema="newick")

    write_fastdate_prior(subtree, dates, var_mult = 0.1, outputfile='node_prior.txt')

    for i in range(int(reps)):
        write_fastdate_tree("unresolved.tre", br_len = 0.01, polytomy_br = 0.001, outputfile='fastdate_input{}.tre'.format(i))
        os.system("fastdate --method_nodeprior --tree_file fastdate_input{}.tre --prior_file node_prior.txt --out_file node_prior{}.tre --max_age {} --bd_rho 1 --grid {}".format(i, i, max_age, max_age*2))

    os.system("sumtrees.py --set-edges=mean-age --summarize-node-ages node_prior*.tre > {}".format(summary))

