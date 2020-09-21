from chronosynth import chronogram



res_ot1000 = map_conflict_ages('ot_1000@tree1')
res_ot1056 = map_conflict_ages('ot_1056@Tr66755')



for synth_node in res_ot1056['supported_nodes']:
    if synth_node in res_ot1000['supported_nodes']:
        print((res_ot1056['supported_nodes'][synth_node]['age'], res_ot1000['supported_nodes'][synth_node]['age']))

# TODO: for test_conf_map_all:
errors = []
no_conflict = []
studies = chronogram.find()
for tag in studies:
    try:
        res = map_conflict_ages(tag)
    except ValueError:
        errors.append(tag)
        print("study", tag, "threw an ERROR on map_conflict ages")
    if res==None:
        print("study", tag, "has no conflict data")
        no_conflict.append(tag)
    else:
        print("study", tag, "has", len(res["supported_nodes"]), "supported_nodes")


studies.index('ot_1876@tree1')
studies[280]
source_id=studies[69]
map_conflict_ages('ot_1876@tree1')
