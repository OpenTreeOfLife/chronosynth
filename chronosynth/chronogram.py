"""Functions for working with chronograms from Phylesystem"""
#!/usr/bin/env python3

import opentree
from opentree import OT

DC = opentree.object_conversion.DendropyConvert()


def find_chronograms():
    """
    Get all study ids in Phylesystem where 'ot:branchLengthMode' == 'ot:time'

    Returns
    -------
    A list of Phylesystem study ids containing at least one tree that is a chronogram.
    """
    output = OT.find_trees(search_property="ot:branchLengthMode", value="ot:time")
    chronograms = set()
    for study in output.response_dict["matched_studies"]:
        study_id = study['ot:studyId']
        for tree in study['matched_trees']:
            tree_id = tree['ot:treeId']
            chronograms.add('{}@{}'.format(study_id, tree_id))
    return list(chronograms)


def get_node_ages(source_id):
    """
    Get node ages for any given tree in the Open Tree of Life

    Returns
    -------
    A list of node ages.
    """
    assert '@' in source_id
    study_id = source_id.split('@')[0]
    tree_id = source_id.split('@')[1]
    study = OT.get_study(study_id)
    study_nexson = study.response_dict['data']
    tree_obj = DC.tree_from_nexson(study_nexson, tree_id)
    time_unit = tree_obj.annotations.get_value("branchLengthTimeUnit")
    assert tree_obj.annotations.get_value("branchLengthMode") == 'ot:time'
    return (tree_obj.internal_node_ages(), time_unit)
