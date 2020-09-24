import chronosynth
from chronosynth.chronogram import find, node_ages, as_dendropy, map_conflict_ages


def test_find():
    studies = find()
    assert 'ot_1000@tree1' in set(studies)
    for tag in studies:
        print(tag)

def test_dp_convert():
    dp_tree = as_dendropy('ot_1000@tree1')


def test_node_ages():
    age_data = node_ages('ot_1000@tree1')
    assert age_data['metadata']['time_unit'] == 'Myr'
    assert len(age_data['ages']) == 290
    assert age_data['ages']['node564'] == 21.381546999999998

def test_conf_map():
    res = map_conflict_ages('ot_1000@tree1')
    print(len(res['supported_nodes']))
    # for i in res['supported_nodes']:
        # print(i)

def test_no_conflict_data():
    res = map_conflict_ages('ot_1980@tree2')

# TODO: add test of map conflict for all from compare_dates.py
def test_conf_map_all():
  pass