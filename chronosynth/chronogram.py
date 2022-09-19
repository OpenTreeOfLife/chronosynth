#!/usr/bin/env python3
"""Functions for working with chronograms from Phylesystem"""

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
import copy

from sh import git
import dendropy

import opentree
from opentree import OT

import peyotl
from peyotl.phylesystem.git_actions import PhylesystemGitAction

import chronosynth

config = configparser.ConfigParser()
config.read(chronosynth.configfile)


#log = logging.getLogger(__name__)
#log.debug("logging set to debug")


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
#    log.debug("api_endpoint is %s", format(OT._api_endpoint))

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
    tree_obj = DC.tree_from_nexson(study_nexson, tree_id, label_format="ot:ottId")
    for tip in tree_obj.leaf_node_iter():
        tip.taxon.label = str(tip.taxon.label)
    return tree_obj

def node_ages(dp_tree,
              ultrametricity_precision=None):
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
    assert dp_tree.annotations.get_value("branchLengthMode") == 'ot:time'
    time_unit = dp_tree.annotations.get_value("branchLengthTimeUnit")
    ages = {}
    try:
        dp_tree.internal_node_ages(ultrametricity_precision=ultrametricity_precision)
        for node in dp_tree.internal_nodes():
            assert node.label not in ages.keys()
            ages[node.label] = node.age
    except dendropy.utility.error.UltrametricityError as err:
        sys.stderr.write("source does not meet ultrametricity_precision threshold of {}".format(ultrametricity_precision))
        return None
    ret = {'ages':ages}
    return ret


def map_conflict_ages(source,
                      compare_to="synth",
                      ultrametricity_precision=None,
                      cache_file_path=None,
                      fresh=False):
    """
    Takes a source id in format study_id@tree_id OR an input tree

    returns a dictionary of:
    {'metadata':{'study_id': study_id, 'tree_id': tree_id,
                 'time_unit': time_unit, 'synth_version':version, sha},
    'supported_nodes':{synth_node_id : {'age':age, 'node_label':node_label}}
    """
    if compare_to == 'synth':
        if cache_file_path is None:
            cache_file_dir = config.get('paths', 'cache_file_dir',
                                        fallback='/tmp/')
            cache_file_path = cache_file_dir + '/{}.json'.format(source)
        if os.path.exists(cache_file_path) and fresh == False:
            sys.stdout.write("Loading" + cache_file_path + "\n")
            ret = json.load(open(cache_file_path))
            return ret
    else:
        cache_file_path = None
    if isinstance(source, str):
        source_id = source
        study_id = source_id.split('@')[0]
        tree_id = source_id.split('@')[1]
        dp_tree = as_dendropy(source_id)
        time_unit = dp_tree.annotations.get_value("branchLengthTimeUnit")
        metadata = {'source_id':source_id,
                    'study_id': study_id,
                    'tree_id': tree_id,
                    'compare_to': compare_to,
                    'time_unit':time_unit}
    elif isinstance(source, dendropy.datamodel.treemodel.Tree):
        dp_tree = source
        metadata = {'source_id':'direct_input',
                    'compare_to': compare_to, 
                    'time_unit':time_unit}
    ages_data = node_ages(dp_tree,
                          ultrametricity_precision=ultrametricity_precision)
    if ages_data is None:
        return None
    version = OT.about()
    metadata['synth_tree_about'] = version['synth_tree_about']

    treestr = conflict_tree_str(dp_tree)
    output_conflict = OT.conflict_str(treestr, compare_to)
    conf = output_conflict.response_dict
    if conf is None:
        sys.stderr.write("No conflict data available for tree {} \n Check its status at {}\n".format(metadata['source_id'],
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
    if cache_file_path:
        ofi = open(cache_file_path, 'w')
        ofi.write(sf)
        ofi.close()
    return ret



def map_conflict_nodes(source, compare_to='synth'):
    """
    Takes a source id in format study_id@tree_id OR a dendropy tree

    returns a dictionary of:
    {
    'matched_nodes':{node_id : synnode_label}}
    }
    """
    if isinstance(source, str):
        source_id = source
        study_id = source_id.split('@')[0]
        tree_id = source_id.split('@')[1]
        dp_tree = as_dendropy(source_id)
        metadata = {'source_id':source_id,'study_id': study_id, 'tree_id': tree_id, 'compare_to': compare_to}
    elif isinstance(source, dendropy.datamodel.treemodel.Tree):
        dp_tree = source
        metadata = {'source_id':'direct_input', 'compare_to': compare_to}
    time_unit = dp_tree.annotations.get_value("branchLengthTimeUnit")
    treestr = conflict_tree_str(dp_tree)
    output_conflict = OT.conflict_str(treestr, compare_to)
    conf = output_conflict.response_dict
    if conf is None:
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


def conflict_tree_str(inputtree):
    """Write out a tree with labels that work for the OpenTree Conflict API
    """
    tmp_tree = copy.deepcopy(inputtree)
    i = 1
    for node in tmp_tree:
        i += 1
        if node.taxon:
            if node.taxon.label.startswith('ott'):
                pass
            else:
                new_label = "ott{}".format(node.taxon.label)
                node.taxon.label = new_label
        else:
            node.label = "{}".format(node.label)
    return tmp_tree.as_string(schema="newick")



def combine_ages_from_sources(source_ids,
                              compare_to = 'synth',
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
            res = map_conflict_ages(tag,
                                    compare_to = compare_to,
                                    ultrametricity_precision=ultrametricity_precision,
                                    fresh=fresh)
        except ValueError:
#            time_unit = res['metadata']['time_unit']
         #   log.info('{}, conflict error\n'.format(tag))
            continue
        if res is None:
 #           log.info('{}, conflict empty\n'.format(tag))
             pass
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


def build_synth_node_source_ages(compare_to='synth',
                                 cache_file_path=None,
                                 ultrametricity_precision=None,
                                 fresh=False,
                                 sources='all'):
    """
    This combines all of the input node ages mapped using "map conflict ages".
    Returns a dictionary, and caches the dict in
    Args:
    cache_file_path (str): Json output. can be given as arg or
                            Defaults to node_ages.json, dir set in config,
                            or defaults /tmp/node_ages.json
    ultrametricity_precision: a float passed to dendropy
f    fresh: Whether to re-map trees to synth and est ages
    """
    if compare_to =='synth':
        if cache_file_path is None:
            cache_file_dir = config.get('paths', 'cache_file_dir',
                                        fallback='/tmp/')
            cache_file_path = cache_file_dir + '/node_ages.json'
        if os.path.exists(cache_file_path) and fresh == False:
            dates = json.load(open(cache_file_path))
            current_sha = get_phylesystem_sha()
            cached_sha = dates['metadata'].get('phylesystem_sha')
            if cached_sha == current_sha:
                sys.stdout.write("No new changes to phylesystem, using cached dates at {}\n".format(cache_file_path))
                return dates
            if cached_sha != current_sha:
                sys.stdout.write("Phylesystem has changed since dates were cached, reloading and saving to {}\n".format(cache_file_path))
        else:
            sys.stdout.write("No date cache found. Loading dates and saving to {}\n".format(cache_file_path))
    else:
        cache_file_path = None
    if sources == 'all':
        sources = find_trees()
    dates = combine_ages_from_sources(sources,
                                      ultrametricity_precision=ultrametricity_precision,
                                      json_out=cache_file_path,
                                      fresh=fresh,
                                      compare_to = compare_to)
    return dates

def synth_node_source_ages(node):
    """
    Return age estimates for a node.
    Arguemnts:
    node: Opentree node id
    """
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
        dates = build_synth_node_source_ages()
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


def date_dist_summary(ages, method, var_mult):
    """Generates a prior for a node based on existing dates"""
    assert method in ['mean', 'max', 'min', 'rand']
    if method == 'mean':
        age_est = sum(ages)/len(ages)
        if len(ages) > 1:
            var = statistics.variance(ages)
        else:
            var = var_mult*age_est
        if var < 0.001:
             var = 0.001 #ToDO haaaaaack
    return(age_est, var)


def write_fastdate_prior(subtree, dates, select='mean', var_mult=0.1, outputfile='node_prior.txt'):
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
                mrca_dict[lab]['tips'] = [str(ti.taxon.label).replace(' ', '_') for ti in nd.leaf_iter()]
    if len(mrca_dict) == 0:
        sys.stderr.write("no calibrations\n")
        return None
    fi = open(outputfile, 'w')
    sources = set()
    for dated_node in mrca_dict:
        fi.write("'")
        fi.write("','".join(mrca_dict[dated_node]['tips']))
        fi.write("'")
        fi.write(' ')
        ages = [source['age'] for source in mrca_dict[dated_node]['ages']]
        sources.update([source['source_id'] for source in mrca_dict[dated_node]['ages']])
        mean, var = date_dist_summary(ages, select, var_mult)
        fi.write('norm({},{},{})\n'.format(0, var, mean))
    fi.close()
    return {'sources':sources, 'outputfile':outputfile}

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


def prune_to_phylo_only(tree, cache_file_dir=None):
    """
    Prune tree to only taxa with some phylogenetic information in OpenTree
    Inputs:
    tree: dendropy formatted tree
    grafted solution: path to current synth grafted solution TODO current default hardcoded hack
    https://files.opentreeoflife.org/synthesis/opentree13.4/output/grafted_solution/grafted_solution.tre
    """
    if cache_file_dir is None:
        cache_file_dir = config.get('paths', 'cache_file_dir',
                                    fallback='/tmp/')
    grafted_solution = cache_file_dir + '/grafted_solution.tre'
    if not os.path.exists(grafted_solution):
        os.system("wget https://files.opentreeoflife.org/synthesis/opentree13.4/output/grafted_solution/grafted_solution.tre -P {}".format(cache_file_dir))
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



def date_synth_subtree(node_id=None,
                       node_ids=None,
                       max_age=None,
                       method = 'bladj',
                       output_dir='Chrono_out',
                       summary='sumtre.tre',
                       phylo_only=False,
                       reps = 5,
                       grid=300,
                       select = 'mean',
                       summarize = False,
                       resolve_polytomies = False):
    """
    Takes a synth subtree subtending from node_id and assigns dates using fastdate.
    Inputs
    ------
    node_id: Opentree synth node id  for root OR
    node_ids: list of ott_ids or node_ids to include
    reps: how many runs to do (polytomies are arbitrarily resolved on each run)
    max_age: maximum age for root node. Default None - will be estimated from data if avail, but is required input if no data
    phylo_only: Prune to only synth tips with phylogenetic information (default False)
    summary: Output. deafult sumtre.tre
    """
    dates = build_synth_node_source_ages(ultrametricity_precision=0.01)
    assert method in ['fastdate','bladj']
    assert select in ['mean', 'random', 'min', 'max']
    assert node_id or node_ids
    if node_ids:
        assert isinstance(node_ids, list)
        root_node = OT.synth_mrca(node_ids=node_ids).response_dict['mrca']['node_id']
    else:
        root_node = node_id

    if max_age:
        max_age_est = max_age
    elif root_node in dates['node_ages']:
        max_age_est = max([source['age'] for source in dates['node_ages'][root_node]]) * 1.25
        sys.stdout.write("Age estimate for root from data: {}\n".format(max_age_est))
    else:
        sys.stderr.write("ERROR: no age estimate for root - please provide max root age using --max_age\n")
        return None
    sys.stdout.write("Root node is {}, age estimate is  {}\n".format(root_node, max_age_est))

    if node_id:
        synth_output = OT.synth_subtree(node_id=root_node, label_format='id', include_all_node_labels=True)
    if node_ids:
        synth_output = OT.synth_induced_tree(node_ids=node_ids, label_format='id', include_all_node_labels=True)
    subtree = dendropy.Tree.get_from_string(synth_output.response_dict['newick'], schema='newick')
    outtreesfile, sources = date_tree(subtree=subtree,
                                      dates=dates,
                                      root_node=root_node,
                                      max_age_est=max_age_est,
                                      method=method,
                                      output_dir=output_dir,
                                      phylo_only=phylo_only,
                                      reps=reps,
                                      grid=grid,
                                      select=select,
                                      resolve_polytomies=resolve_polytomies)
    return_dict = {'dated_trees':open(outtreesfile).readlines(),
                   'topology_sources': synth_output.response_dict['supporting_studies'],
                   'date_sources':sources}
    if summarize:
        summaryfilepath = "{}/{}".format(output_dir, summary)
        sumtree = summarize_trees(outtreesfile, summaryfilepath)
        return_dict['summary_tree':sumtree]
    return return_dict

def date_tree(subtree,
              dates,
              root_node,
              max_age_est,
              method,
              output_dir,
              phylo_only,
              reps,
              resolve_polytomies,
              grid,
              select):
    sys.stdout.write("{} leaves in tree\n".format(len(subtree)))
    if phylo_only:
        subtree = prune_to_phylo_only(subtree)
        sys.stdout.write("{} phylo informed leaves in tree\n".format(len(subtree)))
    if method == 'fastdate':
        outtreesfile, sources = run_fastdate(subtree, dates, max_age_est, output_dir, reps, grid, select)
    if method == 'bladj':
       outtreesfile, sources = run_bladj(subtree, dates, root_node, max_age_est, output_dir, reps=reps, select=select, resolve_polytomies=resolve_polytomies)
    return outtreesfile, sources


def summarize_trees(treesfile, summaryfilepath = "sumtre.tre"):
    treelist = dendropy.TreeList.get(path=treesfile, schema='newick')
    return(summaryfilepath)
     

def run_fastdate(subtree,
                 dates,
                 max_age_est,
                 output_dir,
                 reps,
                 grid,
                 select):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    undated_treefile = "{}/unresolved.tre".format(output_dir)
    subtree.write(path= undated_treefile, schema="newick")
    priorfile = '{}/node_prior.txt'.format(output_dir)
    pr = write_fastdate_prior(subtree, dates, select, var_mult=0.1, outputfile=priorfile)
    if not pr:
        sys.stderr.write("")
    for i in range(int(reps)):
        treefile = '{}/fastdate_input{}.tre'.format(output_dir, i)
        outfile = '{}/fastdate_out{}.tre'.format(output_dir, i)            
        write_fastdate_tree(undated_treefile, br_len=0.01, polytomy_br=0.001, outputfile=treefile)
        fd_cmd = "fastdate --method_nodeprior --tree_file {tf} --prior_file {pf} --out_file {of} --out_form ultrametric --max_age {ma} --bd_rho 1 --grid {gs} > fastdate.out".format(tf=treefile,
                                                                                                                                                               pf=priorfile,
                                                                                                                                                               of=outfile,
                                                                                                                                                               ma=max_age_est,
                                                                                                                                                               gs=grid)
        print(fd_cmd)
        os.system(fd_cmd)
    treesfile = "{}/fastdate_trees.tre".format(output_dir)
    os.system("cat {}/fastdate_out*.tre > {}".format(output_dir, treesfile))
    return treesfile, pr['sources']

def run_bladj(subtree,
              dates,
              root_node,
              max_age_est,
              output_dir,
              select,
              reps,
              resolve_polytomies):
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        undated_treefile = "{}/unresolved.tre".format(output_dir)
        subtree.write(path= undated_treefile, schema="newick")
        for i in range(int(reps)):
            if resolve_polytomies:
                undated_rand_resolve_treefile = "{}/resolved{}.tre".format(output_dir, i)
                subtree.resolve_polytomies()
                subtree.write(path= undated_rand_resolve_treefile, schema="newick")
                pr = write_bladj_ages(subtree, dates, root_node, select, root_age=max_age_est, output_dir=output_dir)
                curr_dir = os.getcwd()
                os.chdir(output_dir)
                os.system("phylocom bladj -f  resolved{}.tre >> bladj.tre".format(i))
                os.chdir(curr_dir)
            else:
                pr = write_bladj_ages(subtree, dates, root_node, select, root_age=max_age_est, output_dir=output_dir)
                curr_dir = os.getcwd()
                os.chdir(output_dir)
                os.system("phylocom bladj -f  unresolved.tre >> bladj.tre")
                os.chdir(curr_dir)
        treesfile = "{}/bladj.tre".format(output_dir)
        return treesfile, pr['sources']

def write_bladj_ages(subtree, dates, root_node, select, root_age, output_dir='.'):
    """
    Writes out an ages and a citation file 
    """
    assert select in ['mean', 'random', 'min', 'max']
    outputfile = "{}/ages".format(output_dir)
    cites = open("{}/date_cites.txt".format(output_dir), 'w')
    ages = open(outputfile,'w')
    dated_nodes = set()
    undated_nodes = set()
    sources = set()
    for node in subtree:
        lab = None
        if node.label:
            if node.label.startswith('mrca'):
                lab = node.label
            elif node.label.startswith('ott'):
                lab = node.label
            if lab in dates['node_ages']:
                dated_nodes.add(lab)
                sources.update([source['source_id'] for source in dates['node_ages'][lab]])
                age_range = [float(source['age']) for source in dates['node_ages'][lab]]
                age_range.sort()
                if select == 'mean':
                    age_est = sum(age_range) / len(age_range) 
                if select == 'min':
                    age_est = age_range[0]
                if select == 'max':
                    age_est = age_range[-1]
                if select == 'random':
                    age_est = random.choice(age_range)
                ages.write("{}\t{}\n".format(node.label, age_est))
            else:
                undated_nodes.add(lab)
    cites.write("\n".join(list(sources)))
    cites.close()
    ages.write("{}\t{}\n".format(root_node, root_age))
    ages.close()
    return {'sources':sources, 'outputfile':outputfile}

def date_custom_synth(custom_synth_dir,
                      method,
                      output_dir,
                      summary,
                      phylo_only,
                      reps,
                      grid):
    dp_tree=dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(custom_synth_dir), schema="newick")
    tips = [leaf.taxon.label for leaf in dp_tree.leaf_node_iter()]
    dates = chronogram.build_synth_node_source_ages(ultrametricity_precision=0.01)
    
    root_node = OT.synth_mrca(node_ids=tips).response_dict['mrca']['node_id']
    if max_age:
        max_age_est = max_age
    elif root_node in dates['node_ages']:
        max_age_est = max([source['age'] for source in dates['node_ages'][root_node]]) * 1.25
        sys.stdout.write("Age estimate for root from data: {}\n".format(max_age_est))
    else:
        sys.stderr.write("ERROR: no age estimate for root - please provide max root age using --max_age\n")
        return None
    sys.stdout.write("Root node is {}, age estimate is  {}\n".format(root_node, max_age_est))


    labeled_tree_str = taxonomy_helpers.conflict_tree_str(dp_tree)
    conf = OT.conflict_str(labeled_tree_str, 'synth').response_dict

    labeled_dp_tree = dendropy.Tree.get_from_string(labeled_tree_str, schema="newick")

    matched_nodes = {}
    for node in labeled_dp_tree:
        if node.label:
            nid = node.label.strip('_')
            if nid in conf:
                node_conf = conf[nid]
                status = node_conf['status']
                witness = node_conf['witness']
                if status == 'supported_by':
                    matched_nodes[nid] = witness
                    node.label = witness
    date_tree()


def get_dated_parent_age(ott_ids=None, node_ids=None, root_node=None):
    '''Tree miust have ott_id + otu labels'''
    dates = build_synth_node_source_ages(ultrametricity_precision=0.01)
    if ott_ids:
        root_node = OT.synth_mrca(ott_ids=ott_ids).response_dict['mrca']['node_id']
    elif node_ids:
        root_node = OT.synth_mrca(node_ids=node_ids).response_dict['mrca']['node_id']
    if root_node is None:
        print("ott_ids, node_ids, or a root_node are required")
        return None     
    if root_node in dates['node_ages']:
        max_age_est = max([source['age'] for source in dates['node_ages'][root_node]])
        return max_age_est
    else:
        ott_id = OT.taxon_mrca(ott_ids=ott_ids).response_dict['mrca']['ott_id']
        if 'ott' + str(ott_id) in dates['node_ages']:
            max_age_est = max([source['age'] for source in dates['node_ages'][root_node]])
            return max_age_est
        else:      
            lineage = OT.taxon_info(ott_id, include_lineage=True).response_dict['lineage']
            for parent in lineage:
                ott_id = parent['ott_id']
                if 'ott' + str(ott_id) in dates['node_ages']:
                    max_age_est = max([source['age'] for source in dates['node_ages']['ott' + str(ott_id)]])
                    print("max_age estimate from {}".format('ott' + str(ott_id)))
                    return max_age_est