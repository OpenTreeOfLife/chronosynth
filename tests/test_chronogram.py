import chronosynth
from chronosynth.chronogram import find_chronograms, get_node_ages


def test_find_chronograms():
    studies = find_chronograms()
    assert 'ot_1000@tree1' in set(studies)

def test_node_ages():
    ages = get_node_ages('ot_1000@tree1')
    assert ages[1] == 'Myr'
    assert len(ages[0]) == 290
