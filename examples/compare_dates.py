from chronosynth import chronogram
from chronosynth.chronogram import find_trees, node_ages, as_dendropy, map_conflict_ages

res_ot1000 = chronogram.map_conflict_ages('ot_1000@tree1')
res_ot1056 = chronogram.map_conflict_ages('ot_1056@Tr66755')



for synth_node in res_ot1056['supported_nodes']:
    if synth_node in res_ot1000['supported_nodes']:
        print((res_ot1056['supported_nodes'][synth_node]['age'], res_ot1000['supported_nodes'][synth_node]['age']))



chronogram.set_dev()


import math

# TODO: for test_conf_map_all:
studies = chronogram.find_trees()

f_errors = open("examples/errors.txt", "w+")
f_no_conflict = open("examples/no_conflict.txt", "w+")
f_good = open("examples/good.txt", "w+")
for tag in studies:
    try:
        res = tag.endswith('tree1')
    except TypeError:
        f_errors.write(tag)
        f_errors.write('\n')
        print("study", tag, "threw an ERROR on mathsqrt")
        continue
    if res==True:
        print("study", tag, "has no conflict data")
        f_no_conflict.write(tag)
        f_no_conflict.write('\n')
    else:
        print("study", tag, "has a different tree number than 1")
        f_good.write(tag)
        f_good.write('\n')

f_errors.close()
f_no_conflict.close()
f_good.close()


f_errors = open("examples/errors.txt", "w+")
f_no_conflict = open("examples/no_conflict.txt", "w+")
f_good = open("examples/good.txt", "w+")
for tag in studies:
    try:
        res = chronogram.map_conflict_ages(tag)
    except ValueError:
        f_errors.write(tag)
        f_errors.write('\n')
        print("study", tag, "threw an ERROR on map_conflict ages")
        continue
    if res==None:
        print("study", tag, "has no conflict data")
        f_no_conflict.write(tag)
        f_no_conflict.write('\n')
    else:
        print("study", tag, "has", len(res["supported_nodes"]), "supported_nodes")
        f_good.write(tag)
        f_good.write('\n')
f_errors.close()
f_no_conflict.close()
f_good.close()

studies.index('ot_1876@tree1')
studies[280]
source_id=studies[69]
map_conflict_ages('ot_1876@tree1')
print((synth_node, res_ot1056['supported_nodes'][synth_node]['age'], res_ot1000['supported_nodes'][synth_node]['age']))
