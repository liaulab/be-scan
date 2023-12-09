"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for generate_BE_guides}
"""

import os
import pytest
from be_scan.sgrna.generate_guides import generate_BE_guides

AR_filepath = 'tests/test_data/sgrna_data/230408_AR_Input.fasta'

@pytest.mark.parametrize("edit_from", ['A', 'C', 'G', 'T'])
@pytest.mark.parametrize("edit_to", ['A', 'C', 'G', 'T'])
@pytest.mark.parametrize("cas_type", ['Sp', 'SpG', 'SpRY'])
@pytest.mark.parametrize("PAM", [None, 'NGG', 'NGN', 'NRN'])
@pytest.mark.parametrize("window", [(4,8), (4,7), (3,9)])
def test_generate_BE_guides_basic_pos(edit_from, edit_to, cas_type, PAM, window):
    df = generate_BE_guides(gene_filepath=AR_filepath,
                            gene_name="AR",
                            edit_from=edit_from, 
                            edit_to=edit_to, 
                            cas_type=cas_type, 
                            window=window,
                            PAM=PAM,
                        )
    assert all(col in df.columns for col in ['sgRNA_seq', 'starting_frame', 'gene_pos', 'chr_pos', 
                                             'exon', 'coding_seq', 'sgRNA_strand', 'gene_strand', 'gene'])
    os.remove('guides.csv')

def test_generate_BE_guides_cas_neg(): 
    # invalid cas_type
    with pytest.raises(Exception): 
        df = generate_BE_guides(gene_filepath=AR_filepath,
                                gene_name="AR",
                                edit_from="C", 
                                edit_to="T", 
                                cas_type='Spp', 
                                window=(4,8),
                                PAM=None,
                        )

def test_generate_BE_guides_window_neg(): 
    # negative window
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(gene_filepath=AR_filepath,
                                gene_name="AR",
                                edit_from="C", 
                                edit_to="T", 
                                cas_type='Sp', 
                                window=(-1, 2),
                                PAM=None,
                        )
    # window backwards
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(gene_filepath=AR_filepath,
                                gene_name="AR",
                                edit_from="C", 
                                edit_to="T", 
                                cas_type='Sp', 
                                window=(4, 1),
                                PAM=None,
                        )

def test_generate_BE_guides_edit_neg(): 
    # edit is not ACGT
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(gene_filepath=AR_filepath,
                                gene_name="AR",
                                edit_from="B", 
                                edit_to="T", 
                                cas_type='Sp', 
                                window=(4,8),
                                PAM=None,
                        )

def test_generate_BE_guides_pam_neg(): 
    # PAM is not valid letters
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(gene_filepath=AR_filepath,
                                gene_name="AR",
                                edit_from="C", 
                                edit_to="T", 
                                cas_type='Sp', 
                                window=(4,8),
                                PAM="ABC",
                        )
        