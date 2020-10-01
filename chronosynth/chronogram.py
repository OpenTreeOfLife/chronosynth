"""Functions for working with chronograms from Phylesystem"""
#!/usr/bin/env python3
import sys
import opentree
import datetime
import json
from opentree import OT
DC = opentree.object_conversion.DendropyConvert()


def set_dev():
    global OT
    OT = opentree.ot_object.OpenTree(api_endpoint='dev')


def set_prod():
    global OT
    OT = opentree.ot_object.OpenTree()

def print_endpoint():
    print(OT._api_endpoint)

def find_trees():
    """
    Get study and tree ids for all chronograms (trees with branch lengths proportional
    to time) in Phylesystem, i.e., where 'ot:branchLengthMode' == 'ot:time'

    Returns
    -------
    A list of Phylesystem chronogram's source ids: study_id@treeid.
    """
    output = OT.find_trees(search_property="ot:branchLengthMode", value="ot:time")
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

def node_ages(source_id):
    """
    Get node ages for a DendroPy tree object.

    Returns
    -------
    A dictionary of dictionaries with
    {'metadata': {'study_id': study_id, 'tree_id': tree_id, 'time_unit': time_unit},
    'ages': {node_label:node_age}
    }

    """
    study_id = source_id.split('@')[0]
    tree_id = source_id.split('@')[1]
    dp_tree = as_dendropy(source_id)
    assert dp_tree.annotations.get_value("branchLengthMode") == 'ot:time'
    time_unit = dp_tree.annotations.get_value("branchLengthTimeUnit")
    metadata = {'study_id': study_id, 'tree_id': tree_id, 'time_unit' :time_unit}
    ages = {}
    dp_tree.internal_node_ages()
    for node in dp_tree.internal_nodes():
        assert node.label not in ages.keys()
        ages[node.label] = node.age
    ret = {'metadata':metadata, 'ages':ages}
    return ret


def map_conflict_ages(source_id):
    """
    Takes a source id in format study_id@tree_id

    returns a dictionary of:
    {'metadata':{'study_id': study_id, 'tree_id': tree_id,
                 'time_unit': time_unit, 'synth_version':version},
    'supported_nodes':{synth_node_id : {'age':age, 'node_label':node_label}}
    }
    """
    ages_data = node_ages(source_id)
    metadata = ages_data['metadata']
    version = OT.about()
    metadata['synth_tree_about'] = version['synth_tree_about']
    output_conflict = OT.conflict_info(study_id=metadata['study_id'], tree_id=metadata['tree_id'])
    conf = output_conflict.response_dict
    if(conf==None):
        url = "https://tree.opentreeoflife.org/curator/study/view/{}/?tab=home&tree={}".format(metadata['study_id'], metadata['tree_id'])
        sys.stderr.write("No conflict data available for tree {} \n Check its status at {}".format(source_id, url))
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


def combine_ages_from_sources(source_ids, json_out = None, failed_sources = 'stderr'):
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
                          'units':'Mya'}
                           ],
                           }

    can optionally write out studies with no conflict response and/or errors to a file speciified by "failed_sources"
    """

    if failed_sources == 'stderr':
        f_errors = sys.stderr
    else: 
        f_errors = open(failed_sources, "w+")

    synth_node_ages = {'metadata':{}, 'node_ages':{}}

    versions = OT.about()
    synth_node_ages['metadata']['synth_tree_about'] = versions['synth_tree_about']
    synth_node_ages['metadata']['date'] = str(datetime.date.today())

    for tag in source_ids:
        try:
            res = map_conflict_ages(tag)
        except ValueError:
            f_errors.write('{}, error\n'.format(tag))
            continue
        if res==None:
            f_errors.write('{}, empty\n'.format(tag))
        else:
            sys.stdout.write("study {} has {} supported nodes\n".format(tag, len(res["supported_nodes"])))
            source_id = tag
            assert synth_node_ages['metadata']['synth_tree_about'] == res['metadata']['synth_tree_about']
            time_unit = res['metadata']['time_unit']
            assert tag == "{}@{}".format(res['metadata']['study_id'],res['metadata']['tree_id'])
            for synth_node in res['supported_nodes']:
                if synth_node not in synth_node_ages['node_ages']: #if not a record yet, we need to create one
                    synth_node_ages['node_ages'][synth_node] = []
                age = res['supported_nodes'][synth_node]['age']
                entry = {'source_id': source_id, 'age':age, 'time_unit':time_unit}
                synth_node_ages['node_ages'][synth_node].append(entry)

    sf = json.dumps(synth_node_ages, sort_keys=True, indent=2, separators=(',', ': '), ensure_ascii=True)
    if json_out is not None:
        ofi = open(json_out,'w')
        ofi.write(sf)
        ofi.close()
    return(synth_node_ages)