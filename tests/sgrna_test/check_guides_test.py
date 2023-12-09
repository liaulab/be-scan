"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for check_guides}
"""

import pytest
from be_scan.sgrna.check_guides import check_guides

genome_file='../reference_genomes/GCF_000001405.26/ncbi_dataset/data/GCF_000001405.26/GCF_000001405.26_GRCh38_genomic.fna'

def test_generate_BE_guides_basic_pos():
    check_guides(guides_file='tests/test_data/sgrna_data/ARSpGCBE_guides.csv', 
                 genome_file=genome_file, 
                 )
    check_guides(guides_file='tests/test_data/sgrna_data/ARSpRYABE_guides.csv', 
                 genome_file=genome_file, 
                 )
    