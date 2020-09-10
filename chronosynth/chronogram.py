"""Functions for working with chronograms from Phylesystem"""
#!/usr/bin/env python3
from opentree import OT


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
        chronograms.add(study_id)
    return list(chronograms)


def get_node_ages(study_id, tree_id):
    """
    Get node ages for any given tree in the Open Tree of Life

    Returns
    -------
    A list of node ages.
    """
    study = OT.get_study(study_id)
    study_nexson = study.response_dict['data']
    tree_obj = DC.tree_from_nexson(study_nexson, tree_id)
    return tree_obj.internal_node_ages()
