"""Functions for working with chronograms from Phylesystem"""
#!/usr/bin/env python3
import sys
import os
import datetime
import json
import subprocess
import re
import random
import configparser
import statistics
import logging

from sh import git
import dendropy

import opentree
from opentree import OT

import peyotl
from peyotl.phylesystem.git_actions import PhylesystemGitAction

import chronosynth

config = configparser.ConfigParser()
config.read(chronosynth.configfile)


log = logging.getLogger(__name__)
log.debug("logging set to debug")


DC = opentree.object_conversion.DendropyConvert()


def set_dev():
    """Set endpoint to dev"""
    global OT
    OT = opentree.ot_object.OpenTree(api_endpoint='dev')


def set_prod():
    """Set endpoint to production"""
    global OT
    OT = opentree.ot_object.OpenTree()

def print_endpoint():
    """Print endpoint"""
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
    Inputs
    ------
    source_id: in format study_id@tree_id
    ultrametricity_precision: passed to dendropy

    Returns
    -------
    A dictionary of dictionaries with
    {'metadata': {'study_id': study_id, 'tree_id': tree_id, 'time_unit': time_unit},
    'ages': {node_label:node_age}
    }

    """
    if ultrametricity_precision is None:
        ultrametricity_precision = float(config.get('params', 'ultrametricity_precision',
                                                    fallback='0.01'))

    study_id = source_id.split('@')[0]
    tree_id = source_id.split('@')[1]
    dp_tree = as_dendropy(source_id)
    assert dp_tree.annotations.get_value("branchLengthMode") == 'ot:time'
    time_unit = dp_tree.annotations.get_value("branchLengthTimeUnit")
    metadata = {'study_id': study_id, 'tree_id': tree_id, 'time_unit' :time_unit}
    ages = {}
    try:
        dp_tree.internal_node_ages(ultrametricity_precision=ultrametricity_precision)
        for node in dp_tree.internal_nodes():
            assert node.label not in ages.keys()
            ages[node.label] = node.age
    except dendropy.utility.error.UltrametricityError as err:
        sys.stderr.write("source {} does not meet ultrametricity_precision threshold of {}".format(source_id, ultrametricity_precision))
        return None
    ret = {'metadata':metadata, 'ages':ages}
    return ret


def map_conflict_ages(source_id,
                      ultrametricity_precision=None,
#                      repo_dir=None, ##ToDo run locally....
                      cache_file_path=None,
                      fresh=False):
    """
    Takes a source id in format study_id@tree_id

    returns a dictionary of:
    {'metadata':{'study_id': study_id, 'tree_id': tree_id,
                 'time_unit': time_unit, 'synth_version':version, sha},
    'supported_nodes':{synth_node_id : {'age':age, 'node_label':node_label}}
    """
    if cache_file_path is None:
        cache_file_dir = config.get('paths', 'cache_file_dir',
                                    fallback='/tmp/')
        cache_file_path = cache_file_dir + '/{}.json'.format(source_id)
    if os.path.exists(cache_file_path) and fresh == False:
        sys.stdout.write("Loading" + cache_file_path + "\n")
        ret = json.load(open(cache_file_path))
        return ret
    ages_data = node_ages(source_id, ultrametricity_precision=ultrametricity_precision)
    if ages_data is None:
        return None
    metadata = ages_data['metadata']
    version = OT.about()
    metadata['synth_tree_about'] = version['synth_tree_about']
    output_conflict = OT.conflict_info(study_id=metadata['study_id'],
                                       tree_id=metadata['tree_id'])
    conf = output_conflict.response_dict
    if conf is None:
        url = "https://tree.opentreeoflife.org/curator/study/view/{}/?tab=home&tree={}".format(metadata['study_id'],
                                                                                               metadata['tree_id'])
        sys.stderr.write("No conflict data available for tree {} \n Check its status at {}\n".format(source_id,
                                                                                                     url))
        return None
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
    sf = json.dumps(ret, sort_keys=True, indent=2, separators=(',', ': '), ensure_ascii=True)
    ofi = open(cache_file_path, 'w')
    ofi.write(sf)
    ofi.close()
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
    assert time_unit == "Mya"
    metadata = {'study_id': study_id, 'tree_id': tree_id}
    output_conflict = OT.conflict_info(study_id, tree_id)
    conf = output_conflict.response_dict
    if conf is None:
        url = "https://tree.opentreeoflife.org/curator/study/view/{}/?tab=home&tree={}".format(metadata['study_id'],
                                                                                               metadata['tree_id'])
        sys.stderr.write("No conflict data available for tree {} \n Check its status at {}\n".format(source_id, url))
        return None
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



def combine_ages_from_sources(source_ids,
                              ultrametricity_precision=None,
                              json_out=None,
                              fresh=False):
    """
    inputs
    ------
    source_ids = a list of source ids in format study_id@tree_id
    ultrametricity_precision = a float passed to dendropy
    json_out = output file
    fresh = if False will re-use cached estimates. If True will re-map all studies.

    Outputs
    -------
    a dictionary
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
            res = map_conflict_ages(tag, ultrametricity_precision=ultrametricity_precision, fresh=fresh)
        except ValueError:
#            time_unit = res['metadata']['time_unit']
            log.info('{}, conflict error\n'.format(tag))
            continue
        if res is None:
            log.info('{}, conflict empty\n'.format(tag))
        else:
            sys.stdout.write("study {} has {} supported nodes\n".format(tag, len(res["supported_nodes"])))
            source_id = tag
            # assert synth_node_ages['metadata']['synth_tree_about'] == res['metadata']['synth_tree_about']
            time_unit = res['metadata']['time_unit']
            if time_unit == 'Myr':
                assert tag == "{}@{}".format(res['metadata']['study_id'], res['metadata']['tree_id']), tag
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
        ofi = open(json_out, 'w')
        ofi.write(sf)
        ofi.close()
    return synth_node_ages


#This should probably go in peyotl or somethings
def get_phylesystem_sha(repo_url="https://github.com/OpenTreeOfLife/phylesystem-1.git", repo_dir=None):
    """Get current phylesystem sha
    """
    if repo_dir:
        assert os.path.exists(repo_dir)
        sha = peyotl.git_storage.git_action.get_HEAD_SHA1('{}/.git'.format(repo_dir))
    else:
        process = subprocess.Popen(["git", "ls-remote", repo_url], stdout=subprocess.PIPE)
        stdout, stderr = process.communicate()
        sha = re.split(r'\t+', stdout.decode('ascii'))[0]
    return sha


def pull_phylesystem(repo_dir, repo_url="https://github.com/OpenTreeOfLife/phylesystem-1.git"):
    """Update local phylesystem"""
    assert os.path.exists(repo_dir)
    git_dir_arg = "--git-dir={}/.git".format(repo_dir)
    git(git_dir_arg, 'pull', repo_url)


def build_synth_node_source_ages(cache_file_path=None, ultrametricity_precision=None, repo_dir=None, fresh=False):
    """
    This combines all of the input node ages mapped using "map conflict ages".
    Returns a dictionary, and caches the dict in
    Args:
    cache_file_path (str): Json output. can be given as arg or
                            Defaults to node_ages.json, dir set in config,
                            or defaults /tmp/node_ages.json
    ultrametricity_precision: a float passed to dendropy
    repo_dir: a local clone of phylesystem. Defaults to None, and uses remote
    fresh: Whether to re-map trees to synth and est ages
    """
    if cache_file_path is None:
        cache_file_dir = config.get('paths', 'cache_file_dir',
                                    fallback='/tmp/')
        cache_file_path = cache_file_dir + '/node_ages.json'
    sources = find_trees()
    if os.path.exists(cache_file_path) and fresh == False:
        dates = json.load(open(cache_file_path))
        current_sha = get_phylesystem_sha(repo_dir=repo_dir)
        cached_sha = dates['metadata'].get('phylesystem_sha')
        if cached_sha == current_sha:
            sys.stdout.write("No new changes to phylesystem, using cached dates at {}\n".format(cache_file_path))
            return dates
        if cached_sha != current_sha:
            if repo_dir:
                repo = PhylesystemGitAction(repo=repo_dir)
                changed_studies = repo.get_changed_docs(cached_sha)
                sys.stdout.write("Mapping changed studies")
                ## Re-est conflict and BL's for changed studies
                changed_trees = [source for source in sources if source.split('@') in changed_studies]
                for source_id in changed_trees:
                    map_conflict_ages(source_id,
                                      ultrametricity_precision=ultrametricity_precision,
                                      repo_dir=repo_dir,
                                      fresh=True)
            else:
                sys.stdout.write("Phylesystem has changed since dates were cached, reloading and saving to {}\n".format(cache_file_path))
    else:
        sys.stdout.write("No date cache found. Loading dates and saving to {}\n".format(cache_file_path))
    dates = combine_ages_from_sources(sources,
                                      ultrametricity_precision=ultrametricity_precision,
                                      json_out=cache_file_path,
                                      fresh=fresh)
    return dates


def synth_node_source_ages(node, cache_file_path=None):
    """
    Return age estimates for a node.
    Arguemnts:
    node: Opentree node id
    cache_file_path: path to a stored json with dates. default None.
    """
    if cache_file_path is None:
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
                retdict['msg'] = msg
        dates = build_synth_node_source_ages(cache_file_path)
        retdict['ot:source_node_ages'] = dates['node_ages'].get(node)
    elif node.startswith('ott') and node.strip('ott').isnumeric():
        tax_resp = OT.taxon_info(node)
        if tax_resp.status_code == 200:
            msg = "Taxon {} is in the taxonomy, but cannot be found in the synth tree\n".format(node)
            retdict['msg'] = msg
            retdict = {'msg': msg, 'synth_response': synth_resp.response_dict, 'tax_response': tax_resp.response_dict}
        else:
            msg = "node {} not found in synthetic tree or taxonomy".format(node)
            retdict = {'msg': msg, 'synth_response': synth_resp.response_dict, 'tax_response': tax_resp.response_dict}
    else:
        msg = "node {} not found in synthetic tree.".format(node)
        retdict = {'msg': msg, 'synth_response': synth_resp.response_dict}
    return retdict


def write_fastdate_prior(subtree, dates, var_mult=0.1, outputfile='node_prior.txt'):
    """
    Writes out a node prior file for fatsdate input, with normal prior on each  node with any dates.
    Where multiple dates exist, variance for normal is estimated from dates.
    Where only one date, variance is date * var_mult.

    Inputs:
    subtree: dendropy tree object labeled with ottids and synth node ids
    dates: dictionary output by synth_node_ages()
    var_mult: Hacky approach to choosing a variance
    outputfile: defaults to node_prior.txt
    """
    mrca_dict = {}
    for node in subtree:
        lab = None
        if node.label:
            lab = str(node.label)
            if lab in dates['node_ages']:
                mrca_dict[lab] = {}
                mrca_dict[lab]['ages'] = dates['node_ages'][lab]
                nd = subtree.find_node_with_label(lab)
                mrca_dict[lab]['tips'] = [ti.taxon.label.replace(' ', '_') for ti in nd.leaf_iter()]
    if len(mrca_dict) == 0:
        sys.stderr.write("no calibrations\n")
        return None
    fi = open(outputfile, 'w')
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

def write_fastdate_tree(subtreepath,
                        br_len=0.01,
                        polytomy_br=0.001,
                        outputfile='fastdate_input.tre'):
    """Takes a subtree from OpenTree synth, with id formatted labels,
    and resolves polytomies and applies arbitrarty branch lengths.
    Uses random to randomies polytomy resolution each time.

    Inputs:
    subtreepath: path to newick tree file
    br_len: branch length to assign to branches
    polytomy_br: branch length to assign to arbitrarily resolved polytomies
    outputfile: default fastdate_input.tre

    """
    subtree = dendropy.Tree.get_from_path(subtreepath, schema="newick")
    subtree.resolve_polytomies(rng=random)
    subtree.suppress_unifurcations()
    for edge in subtree.levelorder_edge_iter():
        if (edge.tail_node is not None) and (edge.length is None):
            edge.length = br_len
        if edge.length == 0:
            edge.length = polytomy_br
    subtree.write(path=outputfile, schema="newick")


def prune_to_phylo_only(tree,
                        grafted_solution="/home/ejmctavish/projects/otapi/opentree13.4_tree/grafted_solution/grafted_solution.tre"):
    """
    Prune tree to only taxa with some phylogenetic information in OpenTree
    Inputs:
    tree: dendropy formatted tree
    grafted solution: path to current synth grafted solution TODO current default hardcoded hack
    https://files.opentreeoflife.org/synthesis/opentree13.4/output/grafted_solution/grafted_solution.tre
    """
    fi = open(grafted_solution).readlines()
    for lin in fi:
        lin = lin.replace('(', ',')
        lin = lin.replace(')', ',')
        lii = lin.split(',')
    synth_ottids = set(lii)
    taxa_to_retain = []
    for leaf in tree.leaf_nodes():
        tax = leaf.taxon
        if tax:
            ottid = tax.label
        else:
            ottid = None
        if ottid in synth_ottids:
            taxa_to_retain.append(tax)
    tree.retain_taxa(taxa_to_retain)
    return tree



def date_synth_subtree(node_id,
                       reps,
                       max_age=None,
                       summary='sumtre.tre',
                       phylo_only=False):
    """
    Takes a synth subtree subtenting from node_id and assigns dates using fastdate.
    Inputs
    ------
    node_id: Opentree synth node id
    reps: how many runs to do (polytomies are arbitrarily resolved on each run)
    max_age: maximum age for root node. Default None - will be estimated from data if avail, but is required input if no data
    phylo_only: Prune to only synth tips with phylogenetic information (default False)
    summary: Output. deafult sumtre.tre
    """
    dates = build_synth_node_source_ages(ultrametricity_precision=0.01)
    if max_age:
        max_age_est = max_age
    elif node_id in dates['node_ages']:
        max_age_est = max([source['age'] for source in dates['node_ages'][node_id]]) * 1.25
    else:
        sys.stderr.write("ERROR: no age estimate for root - please provide max root age using --max_age")
        return None

    output = OT.synth_subtree(node_id=node_id, label_format='id')
    subtree = dendropy.Tree.get_from_string(output.response_dict['newick'], schema='newick')
    sys.stdout.write("{} leaves in tree\n".format(len(subtree)))
    if phylo_only:
        subtree = prune_to_phylo_only(subtree)
        sys.stdout.write("{} phylo informed leaves in tree\n".format(len(subtree)))
    subtree.write(path="unresolved.tre", schema="newick")

    pr = write_fastdate_prior(subtree, dates, var_mult=0.1, outputfile='node_prior.txt')
    if pr:
        for i in range(int(reps)):
            write_fastdate_tree("unresolved.tre", br_len=0.01, polytomy_br=0.001, outputfile='fastdate_input{}.tre'.format(i))
            os.system("fastdate --method_nodeprior --tree_file fastdate_input{}.tre --prior_file node_prior.txt --out_file node_prior{}.tre --max_age {} --bd_rho 1 --grid {} > fastdate.out".format(i,
                                                                                                                                                                                                     i,
                                                                                                                                                                                                     max_age_est,
                                                                                                                                                                                                     max_age_est*2))

        os.system("sumtrees.py --set-edges=mean-age --summarize-node-ages node_prior*.tre > {}".format(summary))
        return summary
    return None

