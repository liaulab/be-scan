"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of integration tests for testing the workflow of findall_BE}
"""
import os
import pytest

# positive control tests

def test_findall_BE_basic_integration_pos(): 
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta CBE SpG test_data/P10275.fasta")
    os.system("rm _SpG_CBE_library.csv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta ABE SpRY test_data/P10275.fasta")
    os.system("rm _SpRY_ABE_library.csv")

def test_findall_BE_opt_output_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta ABE SpRY test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../")
    os.system("rm ../230907_AR_SpRY_ABE_library.csv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta CBE SpRY test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../ --output_type tsv")
    os.system("rm ../230907_AR_SpRY_CBE_library.tsv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta ABE Sp test_data/P10275.fasta --output_prefix 230907_AR")
    os.system("rm 230907_AR_Sp_ABE_library.csv")

def test_findall_BE_opt_pam_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta CBE SpG test_data/P10275.fasta --PAM AGG")
    os.system("rm _SpG_CBE_library.csv")

def test_findall_BE_opt_window_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta CBE SpG test_data/P10275.fasta --window 3 7")
    os.system("rm _SpG_CBE_library.csv")

def test_findall_BE_opt_guide_length_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta CBE SpG test_data/P10275.fasta --guide_length 26") 
    os.system("rm _SpG_CBE_library.csv")


# negative control tests
# guide length can't be 25 but can be 26, error only shows up when it is checking amino acid mutations
