from chronosynth import chronogram
from chronosynth.chronogram import find_trees, node_ages, as_dendropy, map_conflict_ages

# This ones work well:
res_ot1000 = chronogram.map_conflict_ages('ot_1000@tree1')
res_ot1056 = chronogram.map_conflict_ages('ot_1056@Tr66755')

# Compare the dates for the same node
for synth_node in res_ot1056['supported_nodes']:
    if synth_node in res_ot1000['supported_nodes']:
        print((res_ot1056['supported_nodes'][synth_node]['age'], res_ot1000['supported_nodes'][synth_node]['age']))
