import json
from chronosynth import chronogram

#sc = chronogram.find_trees(value='ot:substitutionCount')
#un = chronogram.find_trees(value='ot:undefined')
#cc = chronogram.find_trees(value='ot:changesCount')
#ot = chronogram.find_trees(value='ot:other')

source = 'ot_303@Tr62284'

maps = chronogram.map_conflict_nodes(source)


##Example node prior format
#A,B ln(5,10,50)
#C,D ln(10,3,10)



##Next create priors file using fastdate
dates = json.load(open("node_ages.json"))

mrca_dict = {}
for node_label in maps['matched_nodes']:
    if maps['matched_nodes'][node_label] in dates['node_ages'].keys():
        mrca_dict[node_label] = {}
        mrca_dict[node_label]['synth_node'] = maps['matched_nodes'][node_label]
        mrca_dict[node_label]['ages'] = dates['node_ages'][maps['matched_nodes'][node_label]]
        nd = maps['tree'].find_node_with_label(node_label)
        mrca_dict[node_label]['tips'] = [ti.taxon.label for ti in nd.leaf_iter()]
        assert(maps['tree'].mrca(mrca_dict[node_label]['tips']).label == node_label
