import chronosynth
from chronosynth.chronogram import all_chrono, node_ages, as_dendropy


def test_all():
    studies = all_chrono()
    assert 'ot_1000@tree1' in set(studies)



def test_dp_convert():
    dp_tree = as_dendropy('ot_1000@tree1')


def test_node_ages():
    age_data = node_ages('ot_1000@tree1')
    assert age_data['metadata']['time_unit'] == 'Myr'
    assert len(age_data['ages']) == 290
    assert age_data['ages']['node564'] == 21.381546999999998
