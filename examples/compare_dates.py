from chronosynth.chronogram import map_conflict_ages



res_ot1000 = map_conflict_ages('ot_1000@tree1')
res_ot1056 = map_conflict_ages('ot_1056@Tr66755')



for synth_node in res_ot1056['supported_nodes']:
    if synth_node in res_ot1000['supported_nodes']:
        print((res_ot1056['supported_nodes'][synth_node]['age'], res_ot1000['supported_nodes'][synth_node]['age'])) 
 
