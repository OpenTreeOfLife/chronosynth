
#!/usr/bin/env python3
from opentree import OT


def find_chronograms():
    output = OT.find_trees(search_property="ot:branchLengthMode", value="ot:time")
    chronograms = set()
    for study in output.response_dict["matched_studies"]:
        study_id = study['ot:studyId']
        chronograms.add(study_id)
    return list(chronograms)
