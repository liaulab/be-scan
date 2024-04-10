"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for annotate}
"""

import os
import pytest
from be_scan.sgrna.annotate import annotate

guides1 = "tests/test_data/sgrna/AR_ABESpG_library.csv"
guides2 = "tests/test_data/sgrna/AR_CBESpG_library.csv"
gene_filepath = "tests/test_data/sgrna/230408_AR_Input.fasta"
protein_filepath = "tests/test_data/sgrna/P10275.fasta"

@pytest.mark.parametrize("guides", [guides1, guides2])
@pytest.mark.parametrize("edit_from", ["A", "C", "G", "T"])
@pytest.mark.parametrize("edit_to", ["A", "C", "G", "T"])
@pytest.mark.parametrize("window", [(4,8), (3,9)]) ### amino acids being shifted over by 1 for (4, 7)
@pytest.mark.parametrize("protein", [protein_filepath, ""])
def test_annotate_basic_pos(guides, edit_from, edit_to, window, protein): #, window):
    df = annotate(guides_file      = guides,
                         gene_filepath    = gene_filepath,
                         edit_from        = edit_from,
                         edit_to          = edit_to,
                         protein_filepath = protein, 
                         window           = window,
                         )
    prefix = edit_from+"to"+edit_to
    assert all(col in df.columns for col in ["sgRNA_seq", "starting_frame", "gene_pos", "chr_pos", 
                                             "exon", "coding_seq", "sgRNA_strand", "gene_strand", 
                                             "gene", "domain", 
                                             f"{prefix}_editing_window", f"{prefix}_win_overlap", 
                                             f"{prefix}_target_CDS", f"{prefix}_codon_window", 
                                             f"{prefix}_residue_window", f"{prefix}_edit_site", 
                                             f"{prefix}_mutations", f"{prefix}_muttypes", f"{prefix}_muttype", 
                                             ])
    os.remove("annotated.csv")

def test_annotate_colnames_neg():
    with pytest.raises(Exception): 
        df = annotate(guides_file      = guides1,
                             gene_filepath    = gene_filepath,
                             edit_from        = "C",
                             edit_to          = "T",
                             protein_filepath = protein_filepath, 
                             seq_col = "sgRNA_seqs"
                             )

def test_annotate_edit_neg():
    with pytest.raises(Exception): 
        df = annotate(guides_file      = guides1,
                             gene_filepath    = gene_filepath,
                             edit_from        = "B",
                             edit_to          = "T",
                             protein_filepath = protein_filepath, 
                             )
    with pytest.raises(Exception): 
        df = annotate(guides_file      = guides1,
                             gene_filepath    = gene_filepath,
                             edit_from        = "C",
                             edit_to          = "F",
                             protein_filepath = protein_filepath, 
                             )
    