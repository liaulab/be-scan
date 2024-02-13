"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for reference_check}
"""

import os
import pytest
from be_scan.sgrna.reference_check import reference_check

genome_file='tests/test_data/sgrna/hg38_short.fa'

def test_reference_check_basic_pos():
    reference_check(guides_file='tests/test_data/sgrna/ARSpGCBE_guides.csv', 
                 genome_file=genome_file, 
                 )
    reference_check(guides_file='tests/test_data/sgrna/ARSpRYABE_guides.csv', 
                 genome_file=genome_file, 
                 )
    os.remove('filtered.csv')
    
def test_reference_check_delete_pos():
    df1 = reference_check(guides_file='tests/test_data/sgrna/ARSpGCBE_guides.csv', 
                       genome_file=genome_file, 
                       delete=True
                       )
    df2 = reference_check(guides_file='tests/test_data/sgrna/ARSpGCBE_guides.csv', 
                       genome_file=genome_file, 
                       delete=False
                       )
    assert df1.shape[0] <= df2.shape[0]
    os.remove('filtered.csv')
    