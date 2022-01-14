import json
import os
import dendropy
import random
import sys
from chronosynth import chronogram

from opentree import OTCommandLineTool


def main(arg_list, out, list_for_results=None):
    cli = OTCommandLineTool(usage='Look up studies in the "phylesystem" set of phylogenetic studies that are in the Open '
                        'system')
    cli.parser.add_argument("--ott_id", default=None, required=True,
                        help='The ott_id of the root of the tree.')
    cli.parser.add_argument("--reps", default=10, required=False,
                        help='How many times to randomly resolve polytomies.')
    cli.parser.add_argument("--output", default=None, required=False,
                        help='Output file name.')
    cli.parser.add_argument("--max_age", default=None, required=False,
                        help='max root age.')
    cli.parser.add_argument("--verbose", action="store_true", help='include meta-data in response')
    OT, args = cli.parse_cli(arg_list)


    ##Load dates
    dates = chronogram.build_synth_node_source_ages()

    input_ott_id = args.ott_id
    if 'ott'+input_ott_id in dates['node_ages']:
        max_age = max([source['age'] for source in dates['node_ages']['ott'+input_ott_id]]) * 1.25
    else:
        if args.max_age:
            max_age = args.max_age
        else:
            sys.stderr.write("ERROR: no age estimate for root - please provide max root age using --max_age")
            sys.exit()



    output = OT.synth_subtree(ott_id=input_ott_id, label_format='id')
    subtree = dendropy.Tree.get_from_string(output.response_dict['newick'], schema = 'newick')
    subtree.write(path="unresolved.tre", schema="newick")
    subtree = dendropy.Tree.get_from_path("unresolved.tre", schema = "newick")
    subtree.resolve_polytomies()
    
#                assert(subtree.mrca(taxon_labels=mrca_dict[lab]['tips']).label == lab)
#                Beacuse of our various nodes with outdegree 1, this check not only doesn't work, it alse messes up the tree!






    ##Set some sort of uniform prior on the root informed by taxonomy?
    ## and We actually have a good estimate on bd_rho!!

    for i in range(int(args.reps)):
        subtree=dendropy.Tree.get_from_path("unresolved.tre", schema="newick")
        subtree.resolve_polytomies(rng=random)
        mrca_dict = {}
        internal_node_labels = []
        for node in subtree:
            lab = None
            if node.label:
                lab = str(node.label)
                if lab in dates['node_ages']:
                    mrca_dict[lab] = {}
                    mrca_dict[lab]['ages'] = dates['node_ages'][lab]
                    nd = subtree.find_node_with_label(lab)
                    mrca_dict[lab]['tips'] = [ti.taxon.label.replace(' ','_') for ti in nd.leaf_iter()]
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
        subtree.suppress_unifurcations()
        for edge in subtree.level_order_edge_iter():
            if (edge.tail_node is not None) and (edge.length is None):
                edge.length = 0.01
            if edge.length == 0:
                edge.length = 0.001
        subtree.write(path="test{}.tre".format(i), schema="newick")
        os.system("fastdate --method_nodeprior --tree_file test{}.tre --prior_file test_prior.txt --out_file node_prior{}.tre --max_age {} --bd_rho 1 --grid 100".format(i,i, max_age))


    os.system("sumtrees.py --set-edges=mean-age --summarize-node-ages node_prior*.tre > sumtre.tre")


if __name__  == '__main__':
    rc = main(sys.argv[1:], sys.stdout)
    sys.exit(rc)