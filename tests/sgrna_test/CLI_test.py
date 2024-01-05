# """
# Author: Calvin XiaoYang Hu
# Date: 230911

# {Description: suit of integration tests for testing the workflow of findall_BE}
# """
# import os
# import subprocess
# import filecmp
# import pytest

# # positive control tests
# test_dir = 'tests/test_data/sgrna_data/'
# ref_dir = 'tests/ref_data/sgrna_data/'
# DNA = "230408_AR_Input.fasta"
# protein = "P10275.fasta"

# @pytest.mark.parametrize("query, output, ref", 
#                         [
#                         ("-e1 C -e2 T -c Sp", "_Sp_CtoT_library.csv", "AR_CBESp_library.csv"), # Sp C to T
#                         ("-e1 A -e2 G -c Sp", "_Sp_AtoG_library.csv", "AR_ABESp_library.csv"), # Sp A to G
#                         ("-e1 C -e2 T -c SpG", "_SpG_CtoT_library.csv", "AR_CBESpG_library.csv"), # SpG C to T
#                         ("-e1 A -e2 G -c SpG", "_SpG_AtoG_library.csv", "AR_ABESpG_library.csv"), # SpG A to G
#                         ("-e1 C -e2 T -c SpRY", "_SpRY_CtoT_library.csv", "AR_CBESpRY_library.csv"), # SpRY C to T
#                         ("-e1 A -e2 G -c SpRY", "_SpRY_AtoG_library.csv", "AR_ABESpRY_library.csv"), # SpRY A to G
#                         ("-e1 A -e2 G -c SpRY --output_prefix AR --output_dir ../", "../AR_SpRY_AtoG_library.csv", "AR_ABESpRY_library.csv"), # output A to G
#                         ("-e1 C -e2 T -c SpRY --output_prefix AR --output_dir ../ --output_type tsv", "../AR_SpRY_CtoT_library.tsv", "AR_SpRY_CtoT_library.tsv"), # output C to T
#                         ("-e1 C -e2 T -c SpG --PAM AGG", "_SpG_CtoT_library.csv", "AR_CBESpG_PAMAGG_library.csv"), # pam
#                         ("-e1 C -e2 T -c SpG --window 3 7", "_SpG_CtoT_library.csv", "AR_CBESpG_window37_library.csv"), # window
#                         ]
# )

# def test_findall_BE_basic_integration_posSp(query, output, ref): 
#     out = subprocess.run("python3 -m be_scan findall_be -g {0}{1} -p {0}{2} {3}".format(test_dir, DNA, protein, query), shell=True, capture_output=True)
#     assert out.returncode == 0
#     assert filecmp.cmp(output, "{0}{1}".format(ref_dir, ref))
#     os.system("rm {0}".format(output))


# ### add negative control tests, but don't know how since os.system commands aren't caught when it errors other than with reading the results from "pytest -s"
