"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for generate_BE_guides}
"""

import os
import pytest
from be_scan.sgrna.generate_guides import generate_BE_guides

AR_filepath = "tests/test_data/sgrna/230408_AR_Input.fasta"
AR_params = {"gene_filepath":AR_filepath, "gene_name":"AR"}
CT_Sp_params = {"edit_from":"C", "edit_to":"T", "cas_type":"Sp"}

@pytest.mark.parametrize("edit_from", ["A", "C", "G", "T"])
@pytest.mark.parametrize("edit_to", ["A", "C", "G", "T"])
@pytest.mark.parametrize("cas_type", ["Sp", "SpG", "SpRY"])
def test_generate_BE_guides_basic_required_pos(edit_from, edit_to, cas_type):
    df = generate_BE_guides(**AR_params,
                            edit_from=edit_from, 
                            edit_to=edit_to, 
                            cas_type=cas_type, 
                        )
    assert all(col in df.columns for col in ["sgRNA_seq", "starting_frame", "gene_pos", "chr_pos", 
                                             "exon", "coding_seq", "sgRNA_strand", "gene_strand", "gene"])
    os.remove("guides.csv")

@pytest.mark.parametrize("edit_from", ["A", "C"])
@pytest.mark.parametrize("edit_to", ["G", "T"])
@pytest.mark.parametrize("cas_type", ["SpG"])
@pytest.mark.parametrize("PAM", ["NGN", "NNN", "NGG"])
@pytest.mark.parametrize("window", [[3,7], [5,8]])
@pytest.mark.parametrize("exclude_introns", [True, False])
@pytest.mark.parametrize("exclude_nontargeting", [True, False])
@pytest.mark.parametrize("domains", [{"IDR": [1, 555], "DNA Binding" : [600, 700]}])
def test_generate_BE_guides_basic_optional_pos(edit_from, edit_to, cas_type, 
                                               PAM, window, exclude_introns, exclude_nontargeting, domains):
    df = generate_BE_guides(**AR_params, 
                            edit_from=edit_from, edit_to=edit_to, cas_type=cas_type, 
                            PAM=PAM, 
                            window=window, 
                            exclude_introns=exclude_introns, 
                            exclude_nontargeting=exclude_nontargeting, 
                            domains=domains, 
    )
    assert all(col in df.columns for col in ["sgRNA_seq", "starting_frame", "gene_pos", "chr_pos", 
                                             "exon", "coding_seq", "sgRNA_strand", "gene_strand", 
                                             "gene", "domain"])
    os.remove("guides.csv")

def test_generate_BE_guides_cas_neg(): 
    # invalid cas_type
    with pytest.raises(Exception): 
        df = generate_BE_guides(**AR_params,
                                edit_from="C", edit_to="T", 
                                cas_type="Spp", 
                                window=(4,8), save_df=False
                        )
        
def test_generate_BE_guides_window_neg(): 
    # negative window
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params, **CT_Sp_params, 
                                window=(-1, 2), save_df=False
                        )
    # window backwards
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params, **CT_Sp_params, 
                                window=(4, 1), save_df=False
                        )

def test_generate_BE_guides_edit_neg(): 
    # edit is not ACGT
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params,
                                edit_from="B", edit_to="T", 
                                cas_type="Sp", window=(4,8), save_df=False
                        )
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params,
                                edit_from="C", edit_to="F", 
                                cas_type="Sp", window=(4,8), save_df=False
                        )

def test_generate_BE_guides_pam_neg(): 
    # PAM is not valid letters
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                                PAM="ABC", save_df=False
                        )

def test_generate_BE_guides_excludeintrons_pos(): 
    # make sure that if you exclude introns, dataframe is shorter than when you include them
    df1 = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                             exclude_introns=True, save_df=False
                    )
    df2 = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                             exclude_introns=False, save_df=False
                    )
    assert df1.shape[0] <= df2.shape[0]

def test_generate_BE_guides_excludenontargeting_neg(): 
    # make sure that if you exclude nontargeting, dataframe is shorter than when you include them
    df1 = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                             exclude_nontargeting=True, save_df=False
                    )
    df2 = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                             exclude_nontargeting=False, save_df=False
                    )
    assert df1.shape[0] <= df2.shape[0]

def test_generate_BE_guides_domains_neg(): 
    # what if domains input is not a dict
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                                domains=["A", "B"], save_df=False
                        )
    # what if domains is not string:list of ints
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                                domains={1:[1, 10]}, save_df=False
                        )
    with pytest.raises(AssertionError): 
        df = generate_BE_guides(**AR_params, **CT_Sp_params, window=(4,8),
                                domains={"a":"1-100"}, save_df=False
                        )
    