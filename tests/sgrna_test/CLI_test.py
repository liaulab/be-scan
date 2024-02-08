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
test_dir = 'tests/test_data/sgrna_data/'
ref_dir = 'tests/ref_data/sgrna_data/'
DNA = "230408_AR_Input.fasta"
protein = "P10275.fasta"
genome = "tests/ref_data/sgrna_data/hg38_short.fa"


def test_plot_scatterplot():
    out = subprocess.run("python -m be_scan generate_BE_guides -h", shell=True, capture_output=True)
    assert b"usage: be_scan generate_BE_guides" in out.stdout
    assert b"Generates a list of guides based on a gene .fasta file" in out.stdout

def test_plot_boxes():
    out = subprocess.run("python -m be_scan check_guides -h", shell=True, capture_output=True)
    assert b"usage: be_scan check_guides" in out.stdout
    assert b"Annotates a list of guides with a count of how many times" in out.stdout

def test_plot_corr_heatmap():
    out = subprocess.run("python -m be_scan annotate_guides -h", shell=True, capture_output=True)
    assert b"usage: be_scan annotate_guides" in out.stdout
    assert b"Annotates a list of guides in a dataframe with mutational information" in out.stdout

def test_plot_corr_scatterplot():
    out = subprocess.run("python -m be_scan guides -h", shell=True, capture_output=True)
    assert b"usage: be_scan guides" in out.stdout
    assert b"Generates a list of guides based on a gene .fasta file" in out.stdout


@pytest.mark.parametrize("query, output", 
                        [
                        ("Sp C T", "AR_CBESp_library.csv"), # Sp C to T
                        ("Sp A G", "AR_ABESp_library.csv"), # Sp A to G
                        ("SpG C T", "AR_CBESpG_library.csv"), # SpG C to T
                        ("SpG A G", "AR_ABESpG_library.csv"), # SpG A to G
                        ("SpRY C T", "AR_CBESpRY_library.csv"), # SpRY C to T
                        ("SpRY A G", "AR_ABESpRY_library.csv"), # SpRY A to G
                        ("SpG C T --PAM AGG", "AR_CBESpG_PAMAGG_library.csv"), # pam
                        ("SpG C T --window 3 7", "AR_CBESpG_window37_library.csv"), # window
                        ("SpG C T --exclude_introns", 
                            "AR_CBESpG_exclintron_library.csv"), # SpG C to T
                        ("SpG C T --exclude_nontargeting", 
                            "AR_CBESpG_exclnontarg_library.csv"), # SpG C to T
                        ])
def test_generate_BE_guides_integration_pos(query, output): 
    out = subprocess.run("python3 -m be_scan generate_BE_guides {0}{1} {2} --output_name {3}".format(test_dir, DNA, 
                                                                                                     query, output), 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    assert filecmp.cmp(output, "{0}{1}".format(ref_dir, output))
    # os.system("rm {0}".format(output))


@pytest.mark.parametrize("guides_file, query, output", 
                        [
                        ("AR_CBESpG_library.csv", " C T ", "AR_CBESpG_annotated.csv"), # SpG C to T
                        ("AR_ABESpG_library.csv", " A G ", "AR_ABESpG_annotated.csv"), # SpG A to G
                        ("AR_ABESpG_library.csv", " A G --window 3 7", "AR_ABESpG_annotated.csv"), # SpG A to G
                        ("AR_ABESpG_library.csv", " A G --gene_filepath tests/test_data/sgrna_data/230408_AR_Input.fasta", 
                                                  "AR_ABESpG_annotated.csv"), # SpG A to G
                        ("AR_ABESpG_library.csv", " A G --protein_filepath tests/test_data/sgrna_data/P10275.fasta", 
                                                  "AR_ABESpG_annotated.csv"), # SpG A to G
                        ])
def test_annotate_guides_pos(guides_file, query, output): 
    out = subprocess.run("python3 -m be_scan annotate_guides {0}{1} {2} --output_name {3}".format(ref_dir, guides_file, 
                                                                                                  query, output), 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    # assert filecmp.cmp(output, "{0}{1}".format(ref_dir, output))
    os.system("rm {0}".format(output))
    

@pytest.mark.parametrize("guides_file, output", 
                        [
                        ("AR_CBESpG_library.csv", "AR_CBESpG_checked.csv"), # SpG C to T
                        ("AR_ABESpG_library.csv", "AR_ABESpG_checked.csv"), # SpG A to G
                        ])
def test_check_guides_pos(guides_file, output): 
    out = subprocess.run("python3 -m be_scan check_guides {0}{1} {2} --output_name {3}".format(ref_dir, guides_file, 
                                                                                               genome, output), 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    # assert filecmp.cmp(output, "{0}{1}".format(ref_dir, output))
    os.system("rm {0}".format(output))


@pytest.mark.parametrize("query, output", 
                        [
                        (" Sp C T ", "AR_CBESp_library.csv"), # Sp C to T
                        (" Sp A G ", "AR_ABESp_library.csv"), # Sp A to G
                        (" SpG C T ", "AR_CBESpG_library.csv"), # SpG C to T
                        (" SpG A G ", "AR_ABESpG_library.csv"), # SpG A to G
                        (" SpRY C T ", "AR_CBESpRY_library.csv"), # SpRY C to T
                        (" SpRY A G ", "AR_ABESpRY_library.csv"), # SpRY A to G
                        (" SpG C T --gene_name AR ", "AR_CBESpG_genename_library.csv"), # pam
                        (" SpG C T --protein_filepath tests/test_data/sgrna_data/P10275.fasta ", 
                            "AR_CBESpG_genename_library.csv"), # pam
                        (" SpG C T --PAM AGG ", "AR_CBESpG_PAMAGG_library.csv"), # pam
                        (" SpG C T --window 3 7 ", "AR_CBESpG_window37_library.csv"), # window
                        (" SpG C T --exclude_introns ", 
                            "AR_CBESpG_exclintron_library.csv"), # SpG C to T
                        (" SpG C T --exclude_nontargeting ", 
                            "AR_CBESpG_exclnontarg_library.csv"), # SpG C to T
                        ])
def test_guides_pos(query, output): 
    out = subprocess.run("python3 -m be_scan guides {0}{1} {2} {3} --output_name {4}".format(test_dir, DNA, genome, 
                                                                                             query, output), 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    # assert filecmp.cmp(output, "{0}{1}".format(ref_dir, output))
    os.system("rm {0}".format(output))
