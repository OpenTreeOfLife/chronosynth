import json
import os
import dendropy
import random
import sys
import statistics
from chronosynth import chronogram

from opentree import OTCommandLineTool


def main(arg_list, out, list_for_results=None):
    cli = OTCommandLineTool(usage='Look up studies in the "phylesystem" set of phylogenetic studies that are in the Open '
                        'system')
    cli.parser.add_argument("--ott_id", default=None, required=True,
                        help='The ott_id of the root of the tree.')
    cli.parser.add_argument("--reps", default=10, required=False,
                        help='How many times to randomly resolve polytomies.')
    cli.parser.add_argument("--output", default='sumtre.tre', required=False,
                        help='Output file name.')
    cli.parser.add_argument("--max_age", default=None, required=False,
                        help='max root age.')
    cli.parser.add_argument("--verbose", action="store_true", help='include meta-data in response')
    OT, args = cli.parse_cli(arg_list)

    chronogram.date_synth_subtree(args.ott_id, args.reps, max_age=args.max_age, summary=args.output)

if __name__  == '__main__':
    rc = main(sys.argv[1:], sys.stdout)
    sys.exit(rc)