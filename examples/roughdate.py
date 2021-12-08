import json
import os
from chronosynth import chronogram

sc = chronogram.find_trees(value='ot:substitutionCount')
#un = chronogram.find_trees(value='ot:undefined')
#cc = chronogram.find_trees(value='ot:changesCount')
#ot = chronogram.find_trees(value='ot:other')

#source = 'ot_303@Tr62284'

#source = 'pg_1098@tree2158'
#source = sc[10]

source = 'pg_61@tree816'

maps = chronogram.map_conflict_nodes(source)


##Example node prior format
#A,B ln(5,10,50)
#C,D ln(10,3,10)



##Next create priors file using fastdate
dates = json.load(open("node_ages.json"))

#for taxon in maps['tree'].taxon_namespace:
#     taxon._label = 'ott' + str(taxon.ott_id)

maps['tree'].resolve_polytomies()

for edge in maps['tree'].inorder_edge_iter():
    if (edge.tail_node is not None) and (edge.length is None):
       edge.length = 0.001
    if edge.length == 0:
       edge.length = 0.001



maps['tree'].write(path='test.tre', schema='newick')

mrca_dict = {}
for node_label in maps['matched_nodes']:
    if maps['matched_nodes'][node_label] in dates['node_ages'].keys():
        mrca_dict[node_label] = {}
        mrca_dict[node_label]['synth_node'] = maps['matched_nodes'][node_label]
        mrca_dict[node_label]['ages'] = dates['node_ages'][maps['matched_nodes'][node_label]]
        nd = maps['tree'].find_node_with_label(node_label)
        mrca_dict[node_label]['tips'] = [ti.taxon.label.replace(' ','_') for ti in nd.leaf_iter()]
#        assert(maps['tree'].mrca(taxon_labels=mrca_dict[node_label]['tips']).label == node_label)

##Set some sort of uniform prior on the root informed by taxonomy?
##

fi=open('test_prior.txt','w')
for dated_node in mrca_dict:
    fi.write("'")
    fi.write("','".join(mrca_dict[dated_node]['tips']))
    fi.write("'")
    fi.write(' ')
    ages = [source['age'] for source in mrca_dict[dated_node]['ages']]
    avgage = sum(ages)/len(ages)
    fi.write('ln({},{},{})\n'.format(avgage, avgage, avgage))

fi.close()

os.system("fastdate --method_nodeprior --tree_file test.tre --prior_file test_prior.txt --out_file node_prior.tre --bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 1 --show_tree --grid 100 --rate_mean 5 --rate_variance 1")
