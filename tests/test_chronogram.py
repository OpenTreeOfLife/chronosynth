import chronosynth
from chronosynth.chronogram import all, node_ages


def test_all():
    studies = all()
    assert 'ot_1000@tree1' in set(studies)

def test_node_ages():
    ages = node_ages('ot_1000@tree1')
    assert ages[1] == 'Myr'
    assert len(ages[0]) == 290
