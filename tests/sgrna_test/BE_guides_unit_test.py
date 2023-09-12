"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for testing functions in bigscam/bigscam/sgrna/base_editing_guides}
"""
import sys
import pytest

sys.path.append('../bigscam/sgrna') # relative path to helper directory from test directory
from _gene_ import GeneForCRISPR
from _BE_guides_ import identify_guides

# exhaustive testing of identify_guides function
#   ABE, CBE, CGBE, GCBE, GTBE
#   Sp, SpG, SpRY, SpRY_highefficiency, SpRY_lowefficiency
#   window [i, i+n]: [4,8], [4,7], [3,9], [6,6]
def test_identify_guides_basic_pos(): 
    gene_AR = GeneForCRISPR(filepath='../tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)

    for typ in ['ABE', 'CBE', 'CGBE', 'GCBE', 'GTBE']: 
        for mode in ['Sp', 'SpG', 'SpRY', 'SpRY_highefficiency', 'SpRY_lowefficiency']: 
            for wind in [[4,8], [4,7], [3,9], [6,6]]: 
                fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, 
                                                        cas_type=mode, mode=typ, 
                                                        PAM=None, window=wind)

# window out of range
def test_identify_guides_window_range_neg(): 
    gene_AR = GeneForCRISPR(filepath='../tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='SpG', mode='ABE', 
                                                PAM=None, window=[10, 25])
    
# window negative numbers
def test_identify_guides_window_negative_neg(): 
    gene_AR = GeneForCRISPR(filepath='../tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='SpG', mode='ABE', 
                                                PAM=None, window=[-1, 7])

# cas_type invalid
def test_identify_guides_cas_neg(): 
    gene_AR = GeneForCRISPR(filepath='../tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(Exception): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='SpG', mode='GBE', 
                                                PAM=None, window=[4,7])

# mode invalid
def test_identify_guides_basic_pos(): 
    gene_AR = GeneForCRISPR(filepath='../tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(Exception): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='Spp', mode='ABE', 
                                                PAM=None, window=[4,7])

