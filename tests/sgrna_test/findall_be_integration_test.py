"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of integration tests for testing the workflow of findall_BE}
"""
import os
import subprocess

# positive control tests

def test_findall_BE_basic_integration_pos(): 
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm _SpG_CtoT_library.csv")
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c SpRY -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm _SpRY_AtoG_library.csv")

def test_findall_BE_opt_output_integration_pos():
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c SpRY -p tests/test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../", shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm ../230907_AR_SpRY_AtoG_library.csv")
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpRY -p tests/test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../ --output_type tsv", shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm ../230907_AR_SpRY_CtoT_library.tsv")
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c Sp -p tests/test_data/P10275.fasta --output_prefix 230907_AR", shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm 230907_AR_Sp_AtoG_library.csv")

def test_findall_BE_opt_pam_integration_pos():
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta --PAM AGG", shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm _SpG_CtoT_library.csv")

def test_findall_BE_opt_window_integration_pos():
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta --window 3 7", shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm _SpG_CtoT_library.csv")

def test_findall_BE_opt_guide_length_integration_pos():
    out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta --guide_length 26", shell=True, capture_output=True) 
    assert out.returncode == 0
    os.system("rm _SpG_CtoT_library.csv")
    # out = subprocess.run("python3 -m bigscam findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta --guide_length 25", shell=True, capture_output=True) 
    # assert out.returncode == 0
    # os.system("rm _SpG_CBE_library.csv")


# add negative control tests, but don't know how since os.system commands aren't caught when it errors other than with reading the results from "pytest -s"
