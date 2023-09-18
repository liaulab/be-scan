"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of integration tests for testing the workflow of findall_BE}
"""
import os
import pytest
import filecmp

# positive control tests

def test_findall_BE_basic_integration_pos(): 
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be CBE -c Sp -p test_data/P10275.fasta")
    assert filecmp.cmp("_Sp_CBE_library.csv", "ref_data/AR_CBESp_library.csv")
    os.system("rm _Sp_CBE_library.csv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be ABE -c Sp -p test_data/P10275.fasta")
    assert filecmp.cmp("_Sp_ABE_library.csv", "ref_data/AR_ABESp_library.csv")
    os.system("rm _Sp_ABE_library.csv")

    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be CBE -c SpG -p test_data/P10275.fasta")
    assert filecmp.cmp("_SpG_CBE_library.csv", "ref_data/AR_CBESpG_library.csv")
    os.system("rm _SpG_CBE_library.csv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be ABE -c SpG -p test_data/P10275.fasta")
    assert filecmp.cmp("_SpG_ABE_library.csv", "ref_data/AR_ABESpG_library.csv")
    os.system("rm _SpG_ABE_library.csv")

    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be CBE -c SpRY -p test_data/P10275.fasta")
    assert filecmp.cmp("_SpRY_CBE_library.csv", "ref_data/AR_CBESpRY_library.csv")
    os.system("rm _SpRY_CBE_library.csv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be ABE -c SpRY -p test_data/P10275.fasta")
    assert filecmp.cmp("_SpRY_ABE_library.csv", "ref_data/AR_ABESpRY_library.csv")
    os.system("rm _SpRY_ABE_library.csv")

def test_findall_BE_opt_output_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be ABE -c SpRY -p test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../")
    os.system("rm ../230907_AR_SpRY_ABE_library.csv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be CBE -c SpRY -p test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../ --output_type tsv")
    os.system("rm ../230907_AR_SpRY_CBE_library.tsv")
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be ABE -c Sp -p test_data/P10275.fasta --output_prefix 230907_AR")
    os.system("rm 230907_AR_Sp_ABE_library.csv")

def test_findall_BE_opt_pam_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be CBE -c SpG -p test_data/P10275.fasta --PAM AGG")
    assert filecmp.cmp("_SpG_CBE_library.csv", "ref_data/AR_CBESpG_PAMAGG_library.csv")
    os.system("rm _SpG_CBE_library.csv")

def test_findall_BE_opt_window_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be CBE -c SpG -p test_data/P10275.fasta --window 3 7")
    assert filecmp.cmp("_SpG_CBE_library.csv", "ref_data/AR_CBESpG_window37_library.csv")
    os.system("rm _SpG_CBE_library.csv")

def test_findall_BE_opt_guide_length_integration_pos():
    os.system("python3 ../bigscam/sgrna/findall_BE.py -g test_data/230408_AR_Input.fasta -be CBE -c SpG -p test_data/P10275.fasta --guide_length 26") 
    # assert filecmp.cmp("_SpG_CBE_library.csv", "ref_data/AR_CBESpG_len26_library.csv")
    os.system("rm _SpG_CBE_library.csv")
    # os.system("python3 ../bigscam/sgrna/findall_BE.py test_data/230408_AR_Input.fasta CBE SpG test_data/P10275.fasta --guide_length 25") 
    # os.system("rm _SpG_CBE_library.csv")

# change data format of output to have edits like T878A etc
# add negative control tests, but don't know how since os.system commands aren't caught when it errors other than with reading the results from "pytest -s"
# guide length can't be 25 but can be 26, error only shows up when it is checking amino acid mutations so wasn't caught by unit tests
