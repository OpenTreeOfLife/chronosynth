import chronosynth
import dendropy
from opentree import OT
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
    dp_tree=chronogram.as_dendropy('ot_1000@tree1')
    age_data = node_ages(dp_tree)
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
    assert res['supported_nodes'] == {}
    ## all tips in tree map to same taxon

# TODO: add test of map conflict for all from compare_dates.py
def test_conf_map_all():
    sources = ['ot_1000@tree1','ot_1056@Tr66755']
    resp = chronogram.combine_ages_from_sources(sources)
    assert list(resp.keys()) == ['metadata','node_ages']
    assert len(resp['node_ages']['mrcaott129303ott149204']) == 2
    assert list(resp['node_ages']['mrcaott129303ott149204'][0].keys()) == ['source_id', 'age', 'source_node']

def test_get_phylesystem_sha():
    sha = chronogram.get_phylesystem_sha()
    assert len(sha) == 40

def test_synth_node_source_ages():
    # Hmmmmmm should ideally not require rebuild of whole dang thing...
    ## how to test sha check...
    # Normal synth node
    resp1 = chronogram.synth_node_source_ages('mrcaott1000311ott3643729')

    # Good tax ID, non monophyletic
    resp2 = chronogram.synth_node_source_ages('ott372706')

    # Bad tax ID
    resp3 = chronogram.synth_node_source_ages('ott3727069999999')

    # Bad node id
    resp3 = chronogram.synth_node_source_ages('mrcaott1000311ott364372913412341')



def test_conflict_newick():
    custom_synth_tree = dendropy.Tree.get_from_path("tests/data/labelled_supertree.tre",
                                                    schema = "newick")
    custom_str = chronogram.conflict_tree_str(custom_synth_tree)
    dp_tree = as_dendropy('ot_2013@tree8')
    alt_str = chronogram.conflict_tree_str(dp_tree)
    resp = OT.conflict_str(alt_str, compare_to=custom_str)
    assert len(resp.response_dict.keys()) == 351




def test_custom_synth_node_source_ages():
    custom_synth_tree = dendropy.Tree.get_from_path("tests/data/labelled_supertree.tre",
                                                    schema = "newick")
    custom_str = chronogram.conflict_tree_str(custom_synth_tree)
    custom_dates = chronogram.build_synth_node_source_ages(compare_to = custom_str,
                                                          fresh = True,
                                                          sources = ['ot_2013@tree8'])
    assert len(custom_dates['node_ages']) == 131



def test_fastdate_write():
    # Hmmmmmm should ideally not require rebuild of whole dang thing...
    ## how to test sha check...
    # Normal synth node
    pass