#!/usr/bin/env python3
import dendropy
from chronosynth import chronogram
from opentree import OT, taxonomy_helpers

#Get custom synth:

#!curl -X GET https://ot38.opentreeoflife.org/v3/tree_of_life/custom_built_tree/snacktavish_woodpeckers_1020138_tmp378jz32z.tar.gz --output woodpecker_synth.tar.gz 

custom_synth_dir = "/home/ejmctavish/projects/otapi/OpenTreeCLO/custom_synth_runs/multi_snacktavish_woodpeckers_1020138_tmp39mwh_av"


dp_tree=dendropy.Tree.get_from_path("{}/labelled_supertree/labelled_supertree.tre".format(custom_synth_dir), schema="newick")
tips = [leaf.taxon.label for leaf in dp_tree.leaf_node_iter()]
#Estimate dates a la roughdate
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



dates = chronogram.build_synth_node_source_ages(ultrametricity_precision=0.01)

root_node = OT.synth_mrca(node_ids=tips).response_dict['mrca']['node_id']