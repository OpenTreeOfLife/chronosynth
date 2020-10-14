# TEST manually getting conflict data for all studies with branchLengthMode=="time"

from opentree import OT
from chronosynth import chronogram
import sys
import json
# chronogram.set_dev()
# chronogram.set_prod()

## We don't wnat to run these repeatedly...

outfile = "node_ages1.json"
sources = chronogram.find_trees()
resp = chronogram.combine_ages_from_sources(sources, json_out = outfile, failed_sources='no_conf.txt')


old_resp = json.load(open("node_ages.json"))
sys.stdout.write("Date estimates for {} nodes\n written to {}".format(len(resp["node_ages"]), outfile))


sys.stdout.write("Previous file had date estimates for {} nodes\n written to {}".format(len(old_resp["node_ages"])))
