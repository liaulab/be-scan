"""
Author: Calvin XiaoYang Hu
Date: 230908

{Description: suit of unit tests for testing functions in bigscam/bigscam/sgrna/gene focus is the GeneForCRISPR object}
"""

import pytest

from bigscam.sgrna._gene_ import GeneForCRISPR

def GeneForCRISPR_assert_exons(gene_object): 
    assert isinstance(gene_object.file_content, str)
    assert len(gene_object.exons) == len(gene_object.exons_extra)
    for i in range(len(gene_object.exons)): 
        assert len(gene_object.exons[i]) == len(gene_object.exons_extra[i])-40 
        assert all(c in 'acgtACGT' for c in gene_object.exons[i])
        assert all(c in 'acgtACGT' for c in gene_object.exons_extra[i])

def GeneForCRISPR_assert_guides(gene_object): 
    # assert n is int
    assert len(gene_object.fwd_guides) == len(gene_object.rev_guides)
    for g in gene_object.fwd_guides: assert len(g[0]) == gene_object.n

def workflow(filepath, n): 
    gene_AR = GeneForCRISPR(filepath=filepath)
    gene_AR.parse_exons()
    GeneForCRISPR_assert_exons(gene_AR)
    gene_AR.find_all_guides(n=n)
    GeneForCRISPR_assert_guides(gene_AR)

###

def test_GeneForCRISPR_class_base_pos(): # base positive control test
    # test base case
    workflow(filepath='tests/test_data/230408_AR_Input.fasta', n=23)

def test_GeneForCRISPR_class_n_pos(): # change n positive control test
    # test n != 23
    workflow(filepath='tests/test_data/230408_AR_Input.fasta', n=25)

def test_GeneForCRISPR_class_exon_pos(): # 1 exon file positive control test
    # test 1 exon file
    workflow(filepath='tests/test_data/230408_AR_exon1_Input.fasta', n=23)

def test_GeneForCRISPR_class_empty_pos(): # 1 exon file positive control test
    # test empty gene.fasta file
    workflow(filepath='tests/test_data/230408_AR_empty_Input.fasta', n=23)

###

def test_GeneForCRISPR_class_n_neg(): # negative control test
    # test n not an integer
    with pytest.raises(AssertionError): 
        workflow(filepath='tests/test_data/230408_AR_empty_Input.fasta', n='abc')

def test_GeneForCRISPR_class_filepath_neg(): # negative control test
    # test wrong filepath
    with pytest.raises(FileNotFoundError): 
        workflow(filepath='tests/test_data/nonexistent_file.fasta', n='abc')



