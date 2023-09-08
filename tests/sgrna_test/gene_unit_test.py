"""
Author: Calvin XiaoYang Hu
Date: 230908

{Description: suit of unit tests for testing functions in bigscam/bigscam/sgrna/gene focus is the GeneForCRISPR object}
"""

import sys

sys.path.append('../bigscam/sgrna') # relative path to helper directory from test directory
from gene import GeneForCRISPR

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

###

def test_GeneForCRISPR_class_base_pos(): # base positive control test
    # test base case
    gene_AR1 = GeneForCRISPR(filepath='../tests/test_data/230408_AR_Input.fasta')
    gene_AR1.parse_exons()
    GeneForCRISPR_assert_exons(gene_AR1)
    gene_AR1.find_all_guides(n=23)
    GeneForCRISPR_assert_guides(gene_AR1)

def test_GeneForCRISPR_class_n_pos(): # change n positive control test
    # test n != 23
    gene_AR2 = GeneForCRISPR(filepath='../tests/test_data/230408_AR_Input.fasta')
    gene_AR2.parse_exons()
    GeneForCRISPR_assert_exons(gene_AR2)
    gene_AR2.find_all_guides(n=25)
    GeneForCRISPR_assert_guides(gene_AR2)

def test_GeneForCRISPR_class_exon_pos(): # 1 exon file positive control test
    # test 1 exon file
    gene_AR3 = GeneForCRISPR(filepath='../tests/test_data/230408_AR_exon1_Input.fasta')
    gene_AR3.parse_exons()
    GeneForCRISPR_assert_exons(gene_AR3)
    gene_AR3.find_all_guides(n=23)
    GeneForCRISPR_assert_guides(gene_AR3)

def test_GeneForCRISPR_class_empty_pos(): # 1 exon file positive control test
    # test empty gene.fasta file
    gene_AR4 = GeneForCRISPR(filepath='../tests/test_data/230408_AR_empty_Input.fasta')
    gene_AR4.parse_exons()
    GeneForCRISPR_assert_exons(gene_AR4)
    gene_AR4.find_all_guides(n=23)
    GeneForCRISPR_assert_guides(gene_AR4)

###

def test_GeneForCRISPR_class_neg(): # negative control test
    # test n not an integer
    with pytest.raises(AssertionError): 

    # test wrong filepath




