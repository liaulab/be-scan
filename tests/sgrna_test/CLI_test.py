"""
Author: Calvin XiaoYang Hu
Date: 230911

{Description: suit of integration tests for testing the workflow of findall_BE}
"""
import os
import subprocess
import filecmp
import pytest

# positive control tests
test_dir = "tests/test_data/sgrna/"
DNA = "230408_AR_Input.fasta"
protein = "P10275.fasta"
genome = "tests/test_data/sgrna/hg38_short.fa"
protein = "tests/test_data/sgrna/P10275.fasta"
gene = "tests/test_data/sgrna/230408_AR_Input.fasta"
C = "AR_CBE"
A = "AR_ABE"


def test_plot_scatterplot():
    out = subprocess.run("python -m be_scan generate_library -h", shell=True, capture_output=True)
    assert b"usage: be_scan generate_library"in out.stdout
    assert b"Generates a list of guides based on a gene .fasta file"in out.stdout

def test_boxplot():
    out = subprocess.run("python -m be_scan reference_check -h", shell=True, capture_output=True)
    assert b"usage: be_scan reference_check"in out.stdout
    assert b"Annotates a list of guides with a count of how many times"in out.stdout

def test_corr_heatmap():
    out = subprocess.run("python -m be_scan annotate -h", shell=True, capture_output=True)
    assert b"usage: be_scan annotate"in out.stdout
    assert b"Annotates a list of guides in a dataframe with mutational information"in out.stdout

def test_corr_jointplot():
    out = subprocess.run("python -m be_scan design_library -h", shell=True, capture_output=True)
    assert b"usage: be_scan guides"in out.stdout
    assert b"Generates a list of guides based on a gene .fasta file"in out.stdout


suf = "_library.csv"
@pytest.mark.parametrize("query, output", 
                        [
                        ("Sp C T", f"{C}Sp{suf}"), # Sp C to T
                        ("Sp A G", f"{A}Sp{suf}"), # Sp A to G
                        ("SpG C T", f"{C}SpG{suf}"), # SpG C to T
                        ("SpG A G", f"{A}SpG{suf}"), # SpG A to G
                        ("SpRY C T", f"{C}SpRY{suf}"), # SpRY C to T
                        ("SpRY A G", f"{A}SpRY{suf}"), # SpRY A to G
                        ("SpG C T --PAM AGG", f"{C}SpG_PAMAGG{suf}"), # pam
                        ("SpG C T --window 3 7", f"{C}SpG_window37{suf}"), # window
                        ("SpG C T --exclude_introns", f"{C}SpG_exclintron{suf}"), # SpG C to T
                        ("SpG C T --exclude_nontargeting", f"{C}SpG_exclnontarg{suf}"), # SpG C to T
                        ])
def test_generate_library_integration_pos(query, output): 
    out = subprocess.run(f"python3 -m be_scan generate_library {test_dir}{DNA} {query} --output_name {output}", 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp(output, f"{test_dir}{output}")
    os.system("rm {0}".format(output))


@pytest.mark.parametrize("guides_file, query, output", 
                        [
                        (f"{C}SpG{suf}", "C T", f"{C}SpG_annotated.csv"), # SpG C to T
                        (f"{A}SpG{suf}", "A G", f"{A}SpG_annotated.csv"), # SpG A to G
                        (f"{A}SpG{suf}", "A G --window 3 7", f"{A}SpG_annotated.csv"), # SpG A to G
                        (f"{A}SpG{suf}", f"A G --gene_filepath {gene}", f"{A}SpG_annotated.csv"), # SpG A to G
                        (f"{A}SpG{suf}", f"A G --protein_filepath {protein}", f"{A}SpG_annotated.csv"), # SpG A to G
                        ])
def test_annotate_pos(guides_file, query, output): 
    out = subprocess.run(f"python3 -m be_scan annotate {test_dir}{guides_file} {query} --output_name {output}",
                         shell=True, capture_output=True)
    assert out.returncode == 0
    # assert filecmp.cmp(output, "{0}{1}".format(test_dir, output))
    os.system("rm {0}".format(output))
    

@pytest.mark.parametrize("guides_file, output", 
                        [
                        (f"{C}SpG{suf}", f"{C}SpG_checked.csv"), # SpG C to T
                        (f"{A}SpG{suf}", f"{A}SpG_checked.csv"), # SpG A to G
                        ])
def test_reference_check_pos(guides_file, output): 
    out = subprocess.run(f"python3 -m be_scan reference_check {test_dir}{guides_file} {genome} --output_name {output}",
                         shell=True, capture_output=True)
    assert out.returncode == 0
    # assert filecmp.cmp(output, "{0}{1}".format(test_dir, output))
    os.system("rm {0}".format(output))


@pytest.mark.parametrize("query, output", 
                        [
                        ("Sp C T", f"{C}Sp{suf}"), # Sp C to T
                        ("Sp A G", f"{A}Sp{suf}"), # Sp A to G
                        ("SpG C T", f"{C}SpG{suf}"), # SpG C to T
                        ("SpG A G", f"{A}SpG{suf}"), # SpG A to G
                        ("SpRY C T", f"{C}SpRY{suf}"), # SpRY C to T
                        ("SpRY A G", f"{A}SpRY{suf}"), # SpRY A to G
                        ("SpG C T --gene_name AR", f"{C}SpG_genename{suf}"), # pam
                        (f"SpG C T --protein_filepath {protein}", f"{C}SpG_genename{suf}"), # pam
                        ("SpG C T --PAM AGG", f"{C}SpG_PAMAGG{suf}"), # pam
                        ("SpG C T --window 3 7", f"{C}SpG_window37{suf}"), # window
                        ("SpG C T --exclude_introns", f"{C}SpG_exclintron{suf}"), # SpG C to T
                        ("SpG C T --exclude_nontargeting", f"{C}SpG_exclnontarg{suf}"), # SpG C to T
                        ])
def test_guides_pos(query, output): 
    out = subprocess.run(f"python3 -m be_scan design_library {test_dir}{DNA} {genome} {query} --output_name {output}",
                         shell=True, capture_output=True)
    assert out.returncode == 0
    # assert filecmp.cmp(output, "{0}{1}".format(test_dir, output))
    os.system("rm {0}".format(output))
