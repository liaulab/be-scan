"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of integration tests for testing the workflow of findall_BE}
"""
import os
import subprocess
import filecmp

# positive control tests

def test_findall_BE_basic_integration_pos(): 
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c Sp -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_Sp_CtoT_library.csv", "tests/ref_data/AR_CBESp_library.csv")
    os.system("rm _Sp_CtoT_library.csv")
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c Sp -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_Sp_AtoG_library.csv", "tests/ref_data/AR_ABESp_library.csv")
    os.system("rm _Sp_AtoG_library.csv")

    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_SpG_CtoT_library.csv", "tests/ref_data/AR_CBESpG_library.csv")
    os.system("rm _SpG_CtoT_library.csv")
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c SpG -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_SpG_AtoG_library.csv", "tests/ref_data/AR_ABESpG_library.csv")
    os.system("rm _SpG_AtoG_library.csv")

    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpRY -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_SpRY_CtoT_library.csv", "tests/ref_data/AR_CBESpRY_library.csv")
    os.system("rm _SpRY_CtoT_library.csv")
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c SpRY -p tests/test_data/P10275.fasta", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_SpRY_AtoG_library.csv", "tests/ref_data/AR_ABESpRY_library.csv")
    os.system("rm _SpRY_AtoG_library.csv")

def test_findall_BE_opt_output_integration_pos():
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c SpRY -p tests/test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("../230907_AR_SpRY_AtoG_library.csv", "tests/ref_data/AR_ABESpRY_library.csv")
    os.system("rm ../230907_AR_SpRY_AtoG_library.csv")
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpRY -p tests/test_data/P10275.fasta --output_prefix 230907_AR --output_dir ../ --output_type tsv", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("../230907_AR_SpRY_CtoT_library.tsv", "tests/ref_data/230907_AR_SpRY_CtoT_library.tsv")
    os.system("rm ../230907_AR_SpRY_CtoT_library.tsv")
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 A -e2 G -c Sp -p tests/test_data/P10275.fasta --output_prefix 230907_AR", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("230907_AR_Sp_AtoG_library.csv", "tests/ref_data/AR_ABESp_library.csv")
    os.system("rm 230907_AR_Sp_AtoG_library.csv")

def test_findall_BE_opt_pam_integration_pos():
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta --PAM AGG", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_SpG_CtoT_library.csv", "tests/ref_data/AR_CBESpG_PAMAGG_library.csv")
    os.system("rm _SpG_CtoT_library.csv")

def test_findall_BE_opt_window_integration_pos():
    out = subprocess.run("python3 -m be_scan findall_be -g tests/test_data/230408_AR_Input.fasta -e1 C -e2 T -c SpG -p tests/test_data/P10275.fasta --window 3 7", shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp("_SpG_CtoT_library.csv", "tests/ref_data/AR_CBESpG_window37_library.csv")
    os.system("rm _SpG_CtoT_library.csv")

### add negative control tests, but don't know how since os.system commands aren't caught when it errors other than with reading the results from "pytest -s"
