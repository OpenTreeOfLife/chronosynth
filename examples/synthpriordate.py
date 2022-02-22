import json
import os
import dendropy
import random
import sys
import statistics
from chronosynth import chronogram

from opentree import OTCommandLineTool


def main(arg_list, out, list_for_results=None):
    cli = OTCommandLineTool(usage='date a synth tree')
    cli.parser.add_argument("--node_id", default=None, required=False,
                        help='The ott_id or mrca node id of the root of the tree.')
    cli.parser.add_argument("--node_ids",  nargs='+', default=None, required=False,
                        help='The ott_ids or node ids of desired tips in the tree.')
    cli.parser.add_argument("--reps", default=10, required=False,
                        help='How many times to randomly resolve polytomies.')
    cli.parser.add_argument("--output_dir", default='chrono_out', required=False,
                        help='Output file name.')
    cli.parser.add_argument("--max_age", default=None, required=False,
                        help='max root age.')
    cli.parser.add_argument("--phylo_only", default=False, action='store_true', required=False,
                        help='prune to only tips with some phylogenetic information')
    cli.parser.add_argument("--verbose", action="store_true", help='include meta-data in response')
    OT, args = cli.parse_cli(arg_list)


    if args.node_id:
        chronogram.date_synth_subtree(node_id=args.node_id, reps=args.reps, max_age=args.max_age, output_dir=args.output_dir, phylo_only=args.phylo_only)
    if args.node_ids:
        chronogram.date_synth_subtree(node_ids=args.node_ids,reps=args.reps, max_age=args.max_age, output_dir=args.output_dir, phylo_only=args.phylo_only)

if __name__  == '__main__':

    rc = main(sys.argv[1:], sys.stdout)
    sys.exit(rc)