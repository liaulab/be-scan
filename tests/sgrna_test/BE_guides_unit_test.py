"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for testing functions in bigscam/bigscam/sgrna/base_editing_guides}
"""
import pytest

from bigscam.sgrna._gene_ import GeneForCRISPR
from bigscam.sgrna._BE_guides_ import identify_guides

# exhaustive testing of identify_guides function
#   mode: ABE, CBE, CGBE, GCBE, GTBE
#   type: Sp, SpG, SpRY, SpRY_highefficiency, SpRY_lowefficiency
#   window [i, i+n]: [4,8], [4,7], [3,9], [6,6]
#   PAMs: 'NGG', 'NGN', 'NNN', 'NRN', 'NYN'
def test_identify_guides_basic_pos(): 
    gene_AR = GeneForCRISPR(filepath='tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)

    for mode in ['ABE', 'CBE', 'CGBE', 'GCBE', 'GTBE']: 
        for typ in ['Sp', 'SpG', 'SpRY', 'SpRY_highefficiency', 'SpRY_lowefficiency']: 
            for wind in [[4,8], [4,7], [3,9], [6,6]]: 
                for pams in ['NGG', 'NGN', 'NNN', 'NRN', 'NYN']:
                    fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, 
                                                            cas_type=typ, mode=mode, 
                                                            PAM=pams, window=wind)
                    assert len(mod) == 2
                    ### write more assertions

# window out of range
def test_identify_guides_window_range_neg(): 
    gene_AR = GeneForCRISPR(filepath='tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='SpG', mode='ABE', 
                                                PAM=None, window=[10, 25])
    
# window negative numbers
def test_identify_guides_window_negative_neg(): 
    gene_AR = GeneForCRISPR(filepath='tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='SpG', mode='ABE', 
                                                PAM=None, window=[-1, 7])

# cas_type invalid
def test_identify_guides_cas_neg(): 
    gene_AR = GeneForCRISPR(filepath='tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(Exception): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='SpG', mode='GBE', 
                                                PAM=None, window=[4,7])

# mode invalid
def test_identify_guides_basic_neg(): 
    gene_AR = GeneForCRISPR(filepath='tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(Exception): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='Spp', mode='ABE', 
                                                PAM=None, window=[4,7])


# pam invalid
def test_identify_guides_pam_neg(): 
    gene_AR = GeneForCRISPR(filepath='tests/test_data/230408_AR_Input.fasta')
    gene_AR.parse_exons()
    gene_AR.find_all_guides(n=23)
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_guides(gene_object=gene_AR, cas_type='SpG', mode='ABE', 
                                                PAM='ABC', window=[10, 25])
