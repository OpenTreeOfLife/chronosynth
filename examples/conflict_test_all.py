# TEST manually getting conflict data for all studies with branchLengthMode=="time"

from opentree import OT
from chronosynth import chronogram
import sys
import json
import datetime
# chronogram.set_dev()
# chronogram.set_prod()


#sources = chronogram.find_trees()

f_errors = open("examples/conflict_errors.txt", "w+")
f_no_conflict = open("examples/conflict_none.txt", "w+")
f_good = open("examples/coflict_good.txt", "w+")

synth_node_ages = {'metadata':{}, 'node_ages':{}}

#Key: node_id (or ott id) in synth tree
#example
# synth_node_ages['metadata']{'synth_info':12.3, 
#                        'date_processed':Today}
# synth_node_ages['node_ages']:
# {mrcaott123ott456 : {
#                      [
#                      {'source_id':"ot_1000@tree1", 
#                      'age': 48, 
#                      'node_label':node_label
#                      'units':'Mya'}
#                       ],
#                       }


versions = OT.about()
synth_node_ages['metadata']['synth_tree_about'] = versions['synth_tree_about']
synth_node_ages['metadata']['date'] = str(datetime.date.today())

sources = ['ot_1000@tree1','ot_1056@Tr66755']
for tag in sources:
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
        source_id = tag
        assert synth_node_ages['metadata']['synth_tree_about'] == res['metadata']['synth_tree_about']
        time_unit = res['metadata']['time_unit']
        assert tag == "{}@{}".format(res['metadata']['study_id'],res['metadata']['tree_id'])
        for synth_node in res['supported_nodes']:
            if synth_node not in synth_node_ages['node_ages']: #if not a record yet, we need to create one
                synth_node_ages['node_ages'][synth_node] = []
            age = res['supported_nodes'][synth_node]['age']
            entry = {'source_id': source_id, 'age':age, 'time_unit':time_unit}
            synth_node_ages['node_ages'][synth_node].append(entry)



f_errors.close()
f_no_conflict.close()
f_good.close()

sf = json.dumps(synth_node_ages, sort_keys=True, indent=2, separators=(',', ': '), ensure_ascii=True)
sys.stdout.write(sf)