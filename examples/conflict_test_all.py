# TEST manually getting conflict data for all studies with branchLengthMode=="time"

from opentree import OT
from chronosynth import chronogram

# chronogram.set_dev()
# chronogram.set_prod()


sources = chronogram.find_trees()
resp = chronogram.combine_ages_from_sources(sources, json_out = "node_ages.json", failed_sources='no_conf.txt')

