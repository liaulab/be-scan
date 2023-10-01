from bigscam.analysis import count_reads
import pandas as pd
import os
import pytest

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
