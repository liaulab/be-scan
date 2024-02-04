"""
Author: Simon P Shen, Calvin XiaoYang Hu
Date: 231202

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/count_reads}
"""

from be_scan.analysis import count_reads
import pandas as pd
import os
import pytest
import uuid

file_dir = "tests/test_data/analysis_data/"

# baseline
def test_count_reads():
    count_reads(sample_sheet   = file_dir + "guides_sample_sheet.csv", 
                annotated_lib  = file_dir + "guides_ref.csv",
                file_dir       = file_dir,
                out_dir = file_dir,
                )
    df_counts = pd.read_csv(file_dir + "counts_library.csv", index_col=0, header=None).squeeze("columns")
    assert df_counts.loc["AAAAAAAAAAAAAAAAAAAA"] == '1'
    assert df_counts.loc["TTTTTTTTTTTTTTTTTTTT"] == '2'
    assert len(df_counts) == 3
    # clean up
    os.remove(file_dir + "counts.csv")
    os.remove(file_dir + "noncounts.csv")
    os.remove(file_dir + "stats.txt")

# testing different conditions
@pytest.mark.parametrize("query, ref, match", [
    ("G" + "A" * 20, "A" * 20, True),
    ("G" + "T" * 20, "A" * 20, False),
    ("A" * 20, "A" * 20, True),
    ("G" * 20, "G" * 20, True),
    ("G" * 21, "G" * 20, True),
    ])

def test_matching(query, ref, match):
    try:
        # create fake files
        # create fake .fastq file of reads
        test_id = str(uuid.uuid4())
        fname_query = test_id + ".fastq"
        read = "TTGTGGAAAGGACGAAACACC" + query + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGAC"
        with open(file_dir + fname_query, "wt") as fh:
            fh.writelines([
                "@read1\n",
                read + "\n",
                "+\n",
                "F" * len(read) + "\n"])

        # set up reference .csv of guides
        fname_ref = test_id + ".csv"

        with open(file_dir + fname_ref, "wt") as fh:
            fh.writelines(["sgRNA_seq\n", ref + "\n"])
        
        # set up sample sheet of conditions
        sample_sheet = pd.read_csv(file_dir + "guides_sample_sheet.csv")
        sample_sheet.at[0,'fastq_file'] = test_id + ".fastq"
        fname_sample_sheet = test_id + "sample_sheet.csv"
        sample_sheet.to_csv(file_dir + fname_sample_sheet)

        # run count_reads
        count_reads(sample_sheet   = file_dir + fname_sample_sheet, 
                    annotated_lib  = file_dir + fname_ref,
                    file_dir       = file_dir,
                    out_dir = file_dir,
                    )
        # check whether it's counted as a match
        with open(f"{file_dir}counts.csv", "rt") as fh:
            line = fh.readline().rstrip()
            line = fh.readline().rstrip() # read second line
        assert line == f"{ref},{int(match)}"
        
    finally:
        # clean up
        for fname in [file_dir + fname_query,
                      file_dir + fname_ref,
                      file_dir + fname_sample_sheet, 
                      f"{file_dir}counts.csv",
                      f"{file_dir}noncounts.csv",
                      f"{file_dir}stats.txt"]:
            print(fname)
            if os.path.exists(fname):
                os.remove(fname)
