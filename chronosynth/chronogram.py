"""Functions for working with chronograms from Phylesystem"""
#!/usr/bin/env python3

import opentree
from opentree import OT

DC = opentree.object_conversion.DendropyConvert()


def all():
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
    study = OT.get_study(study_id)
    study_nexson = study.response_dict['data']
    tree_obj = DC.tree_from_nexson(study_nexson, tree_id)
    return tree_obj

def node_ages(dendropy_tree):
    """
    Get node ages for a DendroPy tree object.

    Returns
    -------
    A list of node ages and its branch length type: A tuple.

    """
    time_unit = dendropy_tree.annotations.get_value("branchLengthTimeUnit")
    assert dendropy_tree.annotations.get_value("branchLengthMode") == 'ot:time'
    node_labels = []
    for node in dendropy_tree.internal_nodes():
    print(node.label + ", " + str(node.age))
    node_labels.append(node.label)
    return (tree_obj.internal_node_ages(), time_unit)
