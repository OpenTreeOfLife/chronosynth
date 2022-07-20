import json
import os
import dendropy
import random
import sys
import statistics
from chronosynth import chronogram

from opentree import OT, OTCommandLineTool

"""
e.g. python examples/synthpriordate.py --node_ids_file ../GBIF-Biodiverse-OpenTree/data/gymnosperms.csv --output_dir gymno_fastdate --method fastdate --reps 10 --grid 1000  --max_age 385 --phylo_only

"""


def main(arg_list, out, list_for_results=None):
    cli = OTCommandLineTool(usage='date a synth tree')
    cli.parser.add_argument("--node_id", default=None, required=False,
                        help='The ott_id or mrca node id of the root of the tree.')
    cli.parser.add_argument("--node_ids",  nargs='+', default=None, required=False,
                        help='The ott_ids or node ids of desired tips in the tree.')
    cli.parser.add_argument("--node_ids_file", default=None, required=False,
                        help='The ott_ids or node ids of desired tips in the tree, as first or only column in a file, with header "ott_ids"')
    cli.parser.add_argument("--reps", default=10, required=False,
                        help='How many times to randomly resolve polytomies.')
    cli.parser.add_argument("--grid", default=100, required=False,
                            help='How many grid lines.')
    cli.parser.add_argument("--output_dir", default='chrono_out', required=False,
                        help='Output file name.')
    cli.parser.add_argument("--max_age", default=None, required=False,
                        help='max root age.')
    cli.parser.add_argument("--phylo_only", default=False, action='store_true', required=False,
                        help='prune to only tips with some phylogenetic information')
    cli.parser.add_argument("--method", default='fastdate', required=False,
                        help="method to date tree. Currnetly 'fastdate' or 'bladj")
    cli.parser.add_argument("--verbose", action="store_true", help='include meta-data in response')
    OT, args = cli.parse_cli(arg_list)


    if args.node_id:
        chronogram.date_synth_subtree(node_id=args.node_id,
                                      reps=args.reps,
                                      max_age=args.max_age,
                                      output_dir=args.output_dir,
                                      phylo_only=args.phylo_only,
                                      method=args.method)
    elif args.node_ids:
        chronogram.date_synth_subtree(node_ids=args.node_ids,
                                      reps=args.reps,
                                      max_age=args.max_age,
                                      output_dir=args.output_dir,
                                      phylo_only=args.phylo_only,
                                      method=args.method)
    elif args.node_ids_file:
        assert os.path.exists(args.node_ids_file)
        queryfile = open(args.node_ids_file)
        header = queryfile.readline()
        assert header.split(',')[0] == 'ott_id'
        ott_ids=set()
        for lin in queryfile.readlines():
            lii = lin.split(',')
            if lii[0].startswith('ott' or 'mrca'):
                ott_ids.add(lii[0])
            else:
                ott_ids.add('ott'+lii[0])
        ott_ids = list(ott_ids)    
#        print(OT.synth_mrca(node_ids=ott_ids).response_dict['mrca']['node_id'])
        chronogram.date_synth_subtree(node_ids=ott_ids,
                                      reps=args.reps,
                                      max_age=args.max_age,
                                      output_dir=args.output_dir,
                                      phylo_only=args.phylo_only,
                                      grid=args.grid,
                                      method=args.method)

    else:
        sys.stderr("-node_id OR --node_ids OR -node_ids_file are required as an argument")
        




if __name__  == '__main__':

    rc = main(sys.argv[1:], sys.stdout)
    sys.exit(rc)