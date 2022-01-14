import json
import os
import dendropy
import random
from chronosynth import chronogram
from opentree import OT



##Load dates
dates = json.load(open("../node_ages.json"))

input_ott_id = '889366'

assert('ott'+input_ott_id in dates['node_ages'])
max_age = max([source['age'] for source in dates['node_ages'][input_node]]) * 1.25

output = OT.synth_subtree(ott_id='input_ott_id', label_format='id')

#synth_tree = dendropy.Tree.get_from_path('../../opentree13.4_tree/labelled_supertree/labelled_supertree.tre', schema = 'newick')
#tree2 = dendropy.Tree(dp_tree)
#node_of_interest = tree2.find_node_with_label(input_node)
#subtree = dendropy.Tree(seed_node=node_of_interest) 

subtree = dendropy.Tree.get_from_string(output.response_dict['newick'], schema = 'newick')
subtree.write(path="unresolved.tre", schema="newick")

mrca_dict = {}
for node in subtree:
    lab = None
    if node.label:
        if node.label.startswith('mrca'):
            lab = node.label
        elif node.label.startswith('ott'):
            lab = node.label
        else:
            lab = node.label.split()[-1]
        if lab in dates['node_ages']:
            mrca_dict[lab] = {}
            mrca_dict[lab]['ages'] = dates['node_ages'][lab]
            nd = subtree.find_node_with_label(lab)
            mrca_dict[lab]['tips'] = [ti.taxon.label.replace(' ','_') for ti in nd.leaf_iter()]
            assert(subtree.mrca(taxon_labels=mrca_dict[lab]['tips']).label == lab)

##

if len(mrca_dict) == 0:
    print("no calibrations")
    exit()

fi=open('test_prior.txt','w')
for dated_node in mrca_dict:
    fi.write("'")
    fi.write("','".join(mrca_dict[dated_node]['tips']))
    fi.write("'")
    fi.write(' ')
    ages = [source['age'] for source in mrca_dict[dated_node]['ages']]
    avgage = sum(ages)/len(ages)
    fi.write('norm({},{},{})\n'.format(0, 0.01*avgage, 0.8*avgage))

fi.close()


##Set some sort of uniform prior on the root informed by taxonomy?
## and We actually have a good estimate on bd_rho!!

for i in range(10):
    subtree=dendropy.Tree.get_from_path("unresolved.tre", schema="newick")
    subtree.resolve_polytomies(rng=random)
    for edge in subtree.inorder_edge_iter():
        if (edge.tail_node is not None) and (edge.length is None):
            edge.length = 0.01
        if edge.length == 0:
            edge.length = 0.001
    subtree.write(path="test{}.tre".format(i), schema="newick")
    os.system("fastdate --method_nodeprior --tree_file test{}.tre --prior_file test_prior.txt --out_file node_prior{}.tre --max_age {} --bd_rho 1 --grid 100".format(i,i, max_age))


os.system("sumtrees.py --set-edges=mean-age --summarize-node-ages node_prior*.tre > sumtre.tre")



###### Remind self what the calibratrions even meaaaaaannn