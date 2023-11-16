from be_scan.analysis import count_reads
import pandas as pd
import os
import pytest
import uuid

def test_count_reads():
    with pytest.warns(UserWarning, match="missing column"):
        count_reads("tests/test_data/guides_query.fastq", "tests/test_data/guides_ref.csv",
                    out_counts="tests/test_data/counts.csv",
                    out_np="tests/test_data/counts_np.csv",
                    out_stats="tests/test_data/counts_stats.txt")
    df_counts = pd.read_csv("tests/test_data/counts.csv", index_col=0, header=None).squeeze("columns")
    assert df_counts.loc["AAAAAAAAAAAAAAAAAAAA"] == 1
    assert df_counts.loc["TTTTTTTTTTTTTTTTTTTT"] == 2
    assert len(df_counts) == 2
    # clean up
    os.remove("tests/test_data/counts.csv")
    os.remove("tests/test_data/counts_np.csv")
    os.remove("tests/test_data/counts_stats.txt")

@pytest.mark.parametrize("query,ref,match", [
    ("G" + "A" * 20, "A" * 20, True),
    ("G" + "T" * 20, "A" * 20, False),
    ("A" * 20, "A" * 20, True),
    ("G" * 20, "G" * 20, True),
    ("G" * 21, "G" * 20, True),
    ])
def test_matching(query, ref, match):
    try:
        # create fake files
        test_id = str(uuid.uuid4())
        fname_query = os.path.join("tests/test_data", test_id + ".fastq")
        read = "TTGTGGAAAGGACGAAACACC" + query + "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTAAGCTTGGCGTAACTAGATCTTGAGAC"
        with open(fname_query, "wt") as fh:
            fh.writelines([
                "@read1\n",
                read + "\n",
                "+\n",
                "F" * len(read) + "\n"])

        fname_ref = os.path.join("tests/test_data", test_id + ".csv")
        with open(fname_ref, "wt") as fh:
            fh.writelines(["sgRNA_seq\n", ref + "\n"])
        
        # run count_reads
        with pytest.warns(UserWarning, match="missing column"):
            count_reads(fname_query, fname_ref,
                        out_counts=f"tests/test_data/{test_id}counts.csv",
                        out_np=f"tests/test_data/{test_id}counts_np.csv",
                        out_stats=f"tests/test_data/{test_id}counts_stats.txt")
        # check whether it's counted as a match
        with open(f"tests/test_data/{test_id}counts.csv", "rt") as fh:
            line = fh.readline().rstrip()
        assert line == f"{ref},{int(match)}"
    finally:
        # clean up
        for fname in [fname_query,
                      fname_ref,
                      f"tests/test_data/{test_id}counts.csv",
                      f"tests/test_data/{test_id}counts_np.csv",
                      f"tests/test_data/{test_id}counts_stats.txt"]:
            if os.path.exists(fname):
                os.remove(fname)
