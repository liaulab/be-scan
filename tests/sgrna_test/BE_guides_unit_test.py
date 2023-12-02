"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for testing functions in be_scan/be_scan/sgrna/base_editing_guides}
"""
import pytest

from be_scan.sgrna._gene_ import GeneForCRISPR
from be_scan.sgrna._BE_guides_ import identify_BE_guides

gene_AR = GeneForCRISPR(filepath='tests/test_data/sgrna_data/230408_AR_Input.fasta')
gene_AR.parse_exons()
gene_AR.find_all_guides()

# exhaustive testing of identify_BE_guides function for the following coniditons
#   type: Sp, SpG, SpRY, SpRY_highefficiency, SpRY_lowefficiency
#   window [i, i+n]: [4,8], [4,7], [3,9], [6,6]
#   PAMs: 'NGG', 'NGN', 'NNN', 'NRN', 'NYN'
@pytest.mark.parametrize("edit_from", ['A', 'C', 'G', 'T'])
@pytest.mark.parametrize("edit_to", ['A', 'C', 'G', 'T'])
def test_identify_guides_basic_pos(edit_from, edit_to):
    for typ in ['Sp', 'SpG', 'SpRY', 'SpRY_highefficiency', 'SpRY_lowefficiency']: 
        for wind in [[4,8], [4,7], [3,9], [6,6]]: 
            for pams in ['NGG', 'NGN', 'NNN', 'NRN', 'NYN']:
                fwd_res, rev_res, mod = identify_BE_guides(gene_object=gene_AR, 
                                                        cas_type=typ, edit_from=edit_from, edit_to=edit_to, 
                                                        PAM=pams, window=wind)
                assert len(mod) == 2
                print(len(fwd_res), len(rev_res))
                assert len(fwd_res) > 0 and len(rev_res) > 0

# window out of range
def test_identify_guides_window_range_neg(): 
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_BE_guides(gene_object=gene_AR, cas_type='SpG', edit_from='A', edit_to='G', 
                                                PAM=None, window=[10, 25])
    
# window negative numbers
def test_identify_guides_window_negative_neg(): 
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_BE_guides(gene_object=gene_AR, cas_type='SpG', edit_from='A', edit_to='G', 
                                                PAM=None, window=[-1, 7])

# edit mode invalid
def test_identify_guides_cas_neg(): 
    with pytest.raises(Exception): 
        fwd_res, rev_res, mod = identify_BE_guides(gene_object=gene_AR, cas_type='SpG', edit_from='G', edit_to=None, 
                                                PAM=None, window=[4,7])

# cas_type invalid
def test_identify_guides_basic_neg(): 
    with pytest.raises(Exception): 
        fwd_res, rev_res, mod = identify_BE_guides(gene_object=gene_AR, cas_type='Spp', edit_from='A', edit_to='G', 
                                                PAM=None, window=[4,7])

# pam invalid
def test_identify_guides_pam_neg(): 
    with pytest.raises(AssertionError): 
        fwd_res, rev_res, mod = identify_BE_guides(gene_object=gene_AR, cas_type='SpG', edit_from='A', edit_to='G', 
                                                PAM='ABC', window=[10, 25])
