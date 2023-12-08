"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of integration tests for generate_BE_guides}
"""

import os
import pytest
from be_scan.sgrna.guides import guides

AR_filepath = 'tests/test_data/sgrna_data/230408_AR_Input.fasta'
genome_filepath = '../reference_genomes/GCF_000001405.26/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna'
protein_filepath = 'tests/test_data/sgrna_data/P10275.fasta'

# takes about a minute per test, can't do too many tests
@pytest.mark.parametrize("edit_from", ['A', 'C'])
@pytest.mark.parametrize("edit_to", ['G', 'T'])
@pytest.mark.parametrize("cas_type", ['SpG', 'SpRY'])
@pytest.mark.parametrize("window", [(4,8)]) #, (4,7), (3,9)])
@pytest.mark.parametrize("PAM", [None]) #, 'NGG', 'NGN', 'NRN'])
def test_generate_BE_guides_basic_pos(edit_from, edit_to, cas_type, PAM, window):
    df = guides(gene_filepath=AR_filepath,
                genome_file=genome_filepath,
                protein_filepath=protein_filepath,
                gene_name="AR",
                edit_from=edit_from, 
                edit_to=edit_to, 
                cas_type=cas_type, 
                window=window,
                PAM=PAM,
                )
    assert all(col in df.columns for col in ['sgRNA_seq', 'starting_frame', 'gene_pos', 'chr_pos', 
                                             'exon', 'coding_seq', 'sgRNA_strand', 'gene_strand', 'gene', 
                                             'editing_window', 'win_overlap', 
                                             'target_CDS', 'codon_window', 'residue_window', 'edit_site', 
                                             'mutations', 'muttypes', 'muttype', 
                                             ])
    os.remove('annotated_guides.csv')
