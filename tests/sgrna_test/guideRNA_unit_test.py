
import pytest
from be_scan.sgrna._guideRNA_ import filter_guide, filter_repeats, mutation_combos, format_mutation, annotate_mutations, annotate_dual_mutations, categorize_mutations, calc_editing_window, calc_coding_window, calc_target
from be_scan.sgrna._genomic_ import process_PAM

# for C to T, NGN PAM
@pytest.mark.parametrize("guide", [['AGCTAGCTAGCTAGCTAGCTAGA'], 
                                   ['AGCTAGGTAGCTAGCTAGCTAGA'],
                                   ['AGCTAGCTAGCTAGCTAGCTAAA'],
                                   ['AGCTAGGTAGCTAGCTAGCTAAA'],
                                   ])
@pytest.mark.parametrize("PAM_regex", [process_PAM('NGN')])
@pytest.mark.parametrize("PAM", ['NGN'])
@pytest.mark.parametrize("edit", [('C', 'T')])
@pytest.mark.parametrize("window", [(4, 8)])
def test_filter_guide(guide, PAM_regex, PAM, edit, window): 
    res = filter_guide(guide, PAM_regex, PAM, edit, window)
    assert (res == ('G' == guide[0][21] and 'C' in guide[0][3:8]))

def test_filter_repeats(): 
    sample = [['ABA'], ['ABC'], ['ABA']]
    assert len(filter_repeats(sample))

def test_make_mutations():
    mutated = mutation_combos("ATCG", edit=["C", "T"], dir="sense")
    assert set(mutated) == {'ATCG', 'ATTG'}

def test_make_mutations_processive():
    mutated = mutation_combos("ATCCG", edit=["C", "T"], dir="sense")
    assert set(mutated) == {'ATCCG', 'ATTCG', 'ATCTG', 'ATTTG'}

# def test_format_mutation(): 

# def test_annotate_mutations(): 

# def test_annotate_dual_mutations(): 

# def test_categorize_mutations(): 

# def test_calc_editing_window(): 

# def test_calc_coding_window(): 

# def test_calc_target(): 
