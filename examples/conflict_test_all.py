# TEST manually getting conflict data for all studies with branchLengthMode=="time"

from opentree import OT
from chronosynth import chronogram
import json
# chronogram.set_dev()
# chronogram.set_prod()

## We don't wnat to run these repeatedly...
#sources = chronogram.find_trees()
#resp = chronogram.combine_ages_from_sources(sources, json_out = "node_ages.json", failed_sources='no_conf.txt')


resp = json.load("node_ages.json")
