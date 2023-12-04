"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: filtering }
"""

import pandas as pd
import numpy as np
import ctypes

def check_guides(guides_file,
                 genome, 
                 file_dir='',
                ): 
    
    guides = pd.read_csv("guides_file")
    try: 
        assert 'sgRNA_seq' in guides.columns
    except AssertionError: 
        print("Please make sure 'sgRNA_seq' is a column of guides_file")

    np.savetxt(file_dir+"guides.txt", guides.sgRNA_seq)

    lib = ctypes.CDLL('./main.dll')  # or hello.so if on Linux.
    ref_genome_check = lib.ref_genome_check

    ref_genome_check(file_dir+"guides.txt", file_dir+"GRCh38.fna")
    
