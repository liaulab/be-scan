"""
Author: Calvin XiaoYang Hu
Date: 230908

{Description: suit of unit tests for testing functions in be_scan/be_scan/sgrna/gene focus is the GeneForCRISPR object}
"""

import pytest

from be_scan.sgrna._gene_ import GeneForCRISPR

def GeneForCRISPR_assert_exons(gene_object): 
    assert isinstance(gene_object.file_content, str)
    assert len(gene_object.exons) == len(gene_object.exons_extra)
    for i in range(len(gene_object.exons)): 
        assert len(gene_object.exons[i]) == len(gene_object.exons_extra[i])-40 
        assert all(c in "acgtACGT" for c in gene_object.exons[i])
        assert all(c in "acgtACGT" for c in gene_object.exons_extra[i])

def GeneForCRISPR_assert_guides(gene_object): 
    # assert n is int
    assert len(gene_object.fwd_guides) == len(gene_object.rev_guides)
    for g in gene_object.fwd_guides: assert len(g[0]) == gene_object.n-3

def workflow(filepath): 
    gene_AR = GeneForCRISPR(filepath=filepath)
    gene_AR.parse_exons()
    GeneForCRISPR_assert_exons(gene_AR)
    gene_AR.find_all_guides()
    GeneForCRISPR_assert_guides(gene_AR)

###

def test_GeneForCRISPR_class_base_pos(): # base positive control test
    # test base case
    workflow(filepath="tests/test_data/sgrna/230408_AR_Input.fasta")

def test_GeneForCRISPR_class_exon_pos(): # 1 exon file positive control test
    # test 1 exon file
    workflow(filepath="tests/test_data/sgrna/230408_AR_exon1_Input.fasta")

def test_GeneForCRISPR_class_empty_pos(): # 1 exon file positive control test
    # test empty gene.fasta file
    workflow(filepath="tests/test_data/sgrna/230408_AR_empty_Input.fasta")


def test_GeneForCRISPR_class_filepath_neg(): # negative control test
    # test wrong filepath
    with pytest.raises(FileNotFoundError): 
        workflow(filepath="tests/test_data/sgrna/nonexistent_file.fasta")



