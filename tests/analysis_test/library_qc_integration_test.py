"""
Author: Calvin XiaoYang Hu
Date: 230918

{Description: suit of integration tests for testing the workflow of library_qc}
"""

# import os
# import pytest
# import filecmp

# # test for a small library
# def test_library_qc_basic_integration_pos1(): 
#     os.system("python3 bigscam/analysis/library_qc.py -f tests/test_data/CL-XYH-CBE478-2_S218_L001_R1_001.fastq.gz -i tests/test_data/230228_AR_CBE_exon8.csv -o AR_CBE")
#     assert filecmp.cmp("tests/ref_data/AR_CBE_library_count.csv", "AR_CBE_library_count.csv")
#     assert filecmp.cmp("tests/ref_data/AR_CBE_library_statistics.txt", "AR_CBE_library_statistics.txt")
#     assert filecmp.cmp("tests/ref_data/AR_CBE_imperfect_counts.csv", "AR_CBE_imperfect_counts.csv")
#     os.system("rm AR_CBE_library_count.csv")
#     os.system("rm AR_CBE_library_statistics.txt")
#     os.system("rm AR_CBE_imperfect_counts.csv")

# # test for a medium library
# def test_library_qc_basic_integration_pos2(): 
#     os.system("python3 bigscam/analysis/library_qc.py -f tests/test_data/CL-XYH-ABE479-1_S215_L001_R1_001.fastq.gz -i tests/test_data/230228_AR_ABE_full.csv -o AR_ABE")
#     assert filecmp.cmp("tests/ref_data/AR_ABE_library_count.csv", "AR_ABE_library_count.csv")
#     assert filecmp.cmp("tests/ref_data/AR_ABE_library_statistics.txt", "AR_ABE_library_statistics.txt")
#     assert filecmp.cmp("tests/ref_data/AR_ABE_imperfect_counts.csv", "AR_ABE_imperfect_counts.csv")
#     os.system("rm AR_ABE_library_count.csv")
#     os.system("rm AR_ABE_library_statistics.txt")
#     os.system("rm AR_ABE_imperfect_counts.csv")
    
# # using the wrong file so metrics are bad, and comparing the files with the previous case to make sure they are different
# def test_library_qc_basic_integration_neg(): 
#     os.system("python3 bigscam/analysis/library_qc.py -f tests/test_data/CL-XYH-ABE479-2_S217_L001_R1_001.fastq.gz -i tests/test_data/230228_AR_ABE_full.csv -o AR_ABE2")
#     with pytest.raises(AssertionError): 
#         assert filecmp.cmp("tests/ref_data/AR_ABE_library_count.csv", "AR_ABE2_library_count.csv")
#         assert filecmp.cmp("tests/ref_data/AR_ABE_library_statistics.txt", "AR_ABE2_library_statistics.txt")
#         assert filecmp.cmp("tests/ref_data/AR_ABE_imperfect_counts.csv", "AR_ABE2_imperfect_counts.csv")
#     os.system("rm AR_ABE2_library_count.csv")
#     os.system("rm AR_ABE2_library_statistics.txt")
#     os.system("rm AR_ABE2_imperfect_counts.csv")
