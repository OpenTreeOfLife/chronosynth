import chronosynth
from chronosynth.chronogram import find_chronograms


def test_find_chronograms():
    studies = find_chronograms()
    assert 'ot_1100' in set(studies)
