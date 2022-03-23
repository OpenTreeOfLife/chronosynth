import json
import os
from chronosynth import chronogram
from opentree import OT

sc = chronogram.find_trees(value='ot:substitutionCount')
#un = chronogram.find_trees(value='ot:undefined')
#cc = chronogram.find_trees(value='ot:changesCount')
#ot = chronogram.find_trees(value='ot:other')

#source = 'ot_303@Tr62284'

#source = 'pg_1098@tree2158'
source_id = sc[12]

#source = 'pg_61@tree816'

source_id = 'pg_2551@tree6180'

maps = chronogram.map_conflict_nodes(source_id)
cache_file_path = None

if cache_file_path is None:
        cache_file_dir = chronogram.config.get('paths', 'cache_file_dir',
                                    fallback='/tmp/')
        output_dir = cache_file_dir + '/roughdate/{}'.format(source_id)


##Next create priors file using fastdate
dates = chronogram.build_synth_node_source_ages()

for node in maps['tree']:
    if node.label in maps['matched_nodes']:
        node.label = maps['matched_nodes'][node.label]

maps['tree'].resolve_polytomies()

leaves = [leaf.taxon.label for leaf in maps['tree'].leaf_node_iter()]
ott_ids = []
for leaf in leaves:
    try: 
        ott_ids.append(int((leaf)))
    except:
        pass

root_node = OT.synth_mrca(ott_ids=ott_ids).response_dict['mrca']['node_id']

max_age_est = chronogram.get_dated_parent_age(maps['tree'])

chronogram.date_tree(maps['tree'],
                     dates,
                     root_node,
                     max_age_est,
                     method='fastdate',
                     output_dir=output_dir,
                     summary='trial.tre',
                     phylo_only=False,
                     reps=5,
                     grid=len(leaves))




##Set some sort of uniform prior on the root informed by taxonomy?
## and We actually have a good estimate on bd_rho!!

#os.system("fastdate --method_nodeprior --tree_file test.tre --prior_file test_prior.txt --out_file node_prior.tre --max_age 100 --bd_rho 1 --show_tree --grid 100")


###### Remind self what the calibratrions even meaaaaaannn