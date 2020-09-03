"""Functions for working with chronograms from phylesystem"""
#!/usr/bin/env python3
from opentree import OT


def find_chronograms():
    """Gets all study_ids from phylesystem where 'ot:branchLengthMode' == 'ot:time'
    Returns a list of study_ids"""
    output = OT.find_trees(search_property="ot:branchLengthMode", value="ot:time")
    chronograms = set()
    for study in output.response_dict["matched_studies"]:
        study_id = study['ot:studyId']
        chronograms.add(study_id)
    return list(chronograms)
