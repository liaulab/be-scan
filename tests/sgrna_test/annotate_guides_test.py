"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for annotate_guides}
"""

import os
import pytest
from be_scan.sgrna.annotate_guides import annotate_guides

guides1 = 'tests/test_data/sgrna_data/ARSpGCBE_filtered.csv'
guides2 = 'tests/test_data/sgrna_data/ARSpRYABE_filtered.csv'
gene_filepath = 'tests/test_data/sgrna_data/230408_AR_Input.fasta'
protein_filepath = 'tests/test_data/sgrna_data/P10275.fasta'


@pytest.mark.parametrize("guides", [guides1, guides2])
@pytest.mark.parametrize("edit_from", ['A', 'C', 'G', 'T'])
@pytest.mark.parametrize("edit_to", ['A', 'C', 'G', 'T'])
# @pytest.mark.parametrize("window", [(4,8), (4,7), (3,9)]) ### amino acids being shifted over by 1 for (4, 7)
def test_annotate_guides_basic_pos(guides, edit_from, edit_to): #, window):
    df = annotate_guides(guides_file      = guides,
                         gene_filepath    = gene_filepath,
                         edit_from        = edit_from,
                         edit_to          = edit_to,
                         protein_filepath = protein_filepath, 
                        #  window           = window,
                         )
    assert all(col in df.columns for col in ['sgRNA_seq', 'starting_frame', 'gene_pos', 'chr_pos', 
                                             'exon', 'coding_seq', 'sgRNA_strand', 'gene_strand', 'gene', 
                                             'editing_window', 'win_overlap', 
                                             'target_CDS', 'codon_window', 'residue_window', 'edit_site', 
                                             'mutations', 'muttypes', 'muttype', 
                                             ])
    os.remove('annotated.csv')

def test_annotate_guides_colnames_neg():
    with pytest.raises(Exception): 
        df = annotate_guides(guides_file      = guides1,
                             gene_filepath    = gene_filepath,
                             edit_from        = "C",
                             edit_to          = "T",
                             protein_filepath = protein_filepath, 
                             seq_col = 'sgRNA_seqs'
                             )

def test_annotate_guides_edit_neg():
    with pytest.raises(Exception): 
        df = annotate_guides(guides_file      = guides1,
                             gene_filepath    = gene_filepath,
                             edit_from        = "B",
                             edit_to          = "T",
                             protein_filepath = protein_filepath, 
                             )

# def test_annotate_guides_frame_pos():
    