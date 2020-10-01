import chronosynth
from chronosynth import chronogram
from chronosynth.chronogram import find_trees, node_ages, as_dendropy, map_conflict_ages

def test_setup():
    chronogram.set_prod()

def test_find():
    studies = find_trees()
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
## fails on dev bc study doesn't exist
    chronogram.print_endpoint()
    res = map_conflict_ages('ot_1980@tree2')
    assert res is None
    ## all tips in tree map to same taxon

# TODO: add test of map conflict for all from compare_dates.py
def test_conf_map_all():
    sources = ['ot_1000@tree1','ot_1056@Tr66755']
    resp = chronogram.combine_ages_from_sources(sources)
    assert list(resp.keys()) == ['metadata','node_ages']
    assert len(resp['node_ages']['mrcaott129303ott149204']) == 2
    assert resp['node_ages']['mrcaott129303ott149204'][0].keys() == dict_keys(['source_id', 'age', 'time_unit'])
