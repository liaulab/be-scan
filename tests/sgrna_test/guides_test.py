"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of integration tests for generate_BE_guides}
"""

import os
import pytest
from be_scan.sgrna.guides import guides

AR_filepath = 'tests/test_data/sgrna_data/230408_AR_Input.fasta'
genome_filepath = 'tests/ref_data/sgrna_data/hg38_short.fa' # shorter genome file to speed up tests
protein_filepath = 'tests/test_data/sgrna_data/P10275.fasta'

@pytest.mark.parametrize("edit_from", ['A', 'C'])
@pytest.mark.parametrize("edit_to", ['G', 'T'])
@pytest.mark.parametrize("cas_type", ['SpG', 'SpRY'])
@pytest.mark.parametrize("window", [(4,8), (3,9)]) # (4, 7) doesn't work
@pytest.mark.parametrize("PAM", ['NGG', 'NGN', 'NRN'])
@pytest.mark.parametrize("exclude_introns", [True, False])
@pytest.mark.parametrize("exclude_nontargeting", [True, False])
@pytest.mark.parametrize("domains", [{'IDR': [1, 555], 'DNA Binding' : [600, 700]}])
def test_guides_basic_pos(edit_from, edit_to, cas_type, window, PAM,
                          exclude_introns, exclude_nontargeting, domains
                          ):
    df = guides(gene_filepath=AR_filepath,
                genome_file=genome_filepath,
                protein_filepath=protein_filepath,
                gene_name="AR",
                edit_from=edit_from, 
                edit_to=edit_to, 
                cas_type=cas_type, 
                window=window,
                PAM=PAM,
                exclude_introns=exclude_introns, 
                exclude_nontargeting=exclude_nontargeting, 
                domains=domains,
                )
    prefix = edit_from+'to'+edit_to
    assert all(col in df.columns for col in ['sgRNA_seq', 'starting_frame', 'gene_pos', 'chr_pos', 
                                             'exon', 'coding_seq', 'sgRNA_strand', 'gene_strand', 'gene', 
                                             prefix+'_editing_window', prefix+'_win_overlap', 
                                             prefix+'_target_CDS', prefix+'_codon_window', 
                                             prefix+'_residue_window', prefix+'_edit_site', 
                                             prefix+'_mutations', prefix+'_muttypes', prefix+'_muttype', 
                                             ])
    os.remove('annotated_guides.csv')
