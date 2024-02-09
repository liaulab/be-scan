"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of unit tests for check_guides}
"""

import os
import pytest
from be_scan.sgrna.check_guides import check_guides

genome_file='tests/test_data/sgrna/hg38_short.fa'

def test_check_guides_basic_pos():
    check_guides(guides_file='tests/test_data/sgrna/ARSpGCBE_guides.csv', 
                 genome_file=genome_file, 
                 )
    check_guides(guides_file='tests/test_data/sgrna/ARSpRYABE_guides.csv', 
                 genome_file=genome_file, 
                 )
    os.remove('filtered.csv')
    
def test_check_guides_delete_pos():
    df1 = check_guides(guides_file='tests/test_data/sgrna/ARSpGCBE_guides.csv', 
                       genome_file=genome_file, 
                       delete=True
                       )
    df2 = check_guides(guides_file='tests/test_data/sgrna/ARSpGCBE_guides.csv', 
                       genome_file=genome_file, 
                       delete=False
                       )
    assert df1.shape[0] <= df2.shape[0]
    os.remove('filtered.csv')
    