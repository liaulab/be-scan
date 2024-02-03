
import pytest
from be_scan.sgrna._guideRNA_ import *
from be_scan.sgrna._genomic_ import process_PAM

# for C to T, NGN PAM
@pytest.mark.parametrize("guide", [['AGCTAGCTAGCTAGCTAGCT', 'AGA'], 
                                   ['AGCTAGGTAGCTAGCTAGCT', 'AGA'],
                                   ['AGCTAGCTAGCTAGCTAGCT', 'AAA'],
                                   ['AGCTAGGTAGCTAGCTAGCT', 'AAA'],
                                   ])
@pytest.mark.parametrize("PAM_regex", [process_PAM('NGN')])
@pytest.mark.parametrize("edit", [('C', 'T')])
@pytest.mark.parametrize("window", [(4, 8)])
@pytest.mark.parametrize("exclude_introns", [False, True])
@pytest.mark.parametrize("exclude_nontargeting", [False, True])
def test_filter_guide(guide, PAM_regex, edit, window, exclude_introns, exclude_nontargeting): 
    res = filter_guide(guide, PAM_regex, edit, window, exclude_introns, exclude_nontargeting)
    if exclude_nontargeting:
        assert (res == ('G' == guide[1][1] and 'C' in guide[0][3:8]))

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
