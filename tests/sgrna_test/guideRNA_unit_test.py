
import pytest
from be_scan.sgrna._guideRNA_ import *
from be_scan.sgrna._genomic_ import process_PAM
import pandas as pd

# VARIABLES
seq_col = 'sgRNA_seq'
gene_pos_col='gene_pos'
frame_col = 'starting_frame'
strand_col = 'sgRNA_strand'
col_names = frame_col, strand_col, gene_pos_col, seq_col

sample = {'sgRNA_seq':      ['AGCTAGCTAGCTAGCTAGCT', 'AGCTAGCTAGCTAGCTAGCT', 
                             'AGCTAGCTAGCTAGCTAGCT', 'AGCTAGCTAGCTAGCTAGCT', 
                             'AGCTAGCTAGCTAGCTAGCT', 'AGCTAGCTAGCTAGCTAGCT'], 
          'gene_pos':       [10, 10, 10, 10, 10, 10], 
          'starting_frame': [0, 1, 2, 0, 1, 2], 
          'sgRNA_strand':   ['sense', 'antisense', 'sense', 
                             'antisense', 'sense', 'antisense'],
        }
sample_df = pd.DataFrame.from_dict(sample)

# TESTS
# for C to T, NGN PAM
@pytest.mark.parametrize("guide", [['AGCTAGCTAGCTAGCTAGCT', 'AGA'], 
                                   ['AGCTAGGTAGCTAGCTAGCT', 'AGA'],
                                   ['AGCTAGCTAGCTAGCTAGCT', 'AAA'], # no PAM
                                   ['AGCTAGGTAGCTAGCTAGCT', 'AAA'], # no PAM
                                   ['agctaggtAGCTAGCTAGCT', 'AGA'], # no PAM
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

# def test_annotate_mutations(): 

# def test_annotate_dual_mutations(): 

def test_mutation_combos():
    mutated = mutation_combos("ATCG", edit=["C", "T"], dir="sense")
    assert set(mutated) == {'ATCG', 'ATTG'}

def test_mutation_combos_processive():
    mutated = mutation_combos("ATCCG", edit=["C", "T"], dir="sense")
    assert set(mutated) == {'ATCCG', 'ATTCG', 'ATCTG', 'ATTTG'}

def test_format_mutation(): 
    assert format_mutation('MEV', 'MPV', 1, None, "dummy_guide_sequence") == "E2P"

# def test_categorize_mutations(): 
    
# def test_calc_target(): 

# def test_calc_coding_window(): 

@pytest.mark.parametrize("window", [[4,8], [3,7]])
def test_calc_editing_window(window): 
    sample_df['editing_window'] = sample_df.apply(lambda x: calc_editing_window(x, window, col_names), axis=1)

# def test_parse_position_frames(): 

def test_extract_uppercase_letters(): 
    assert extract_uppercase_letters('acgtACGT') == 'ACGT'

def test_find_first_uppercase_index(): 
    assert find_first_uppercase_index('acgtACGT') == 4

def test_annotate_intron_exon(): 
    assert annotate_intron_exon('acgtACGT', [4, 8]) == "Exon/Intron"
    assert annotate_intron_exon('acgtacgt', [4, 8]) == "Intron"
    assert annotate_intron_exon('ACGTACGT', [4, 8]) == "Exon"
