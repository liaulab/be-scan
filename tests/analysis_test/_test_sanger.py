from be_scan.analysis import validate_cloning
from be_scan.analysis._sanger import _golden_gate
import pandas as pd
import pytest
import subprocess
import os

@pytest.mark.parametrize("vector_fname", ["tests/test_data/sanger/pLentiCRISPRv2.gbk", "tests/test_data/sanger/CROPseq-Guide-Puro.gbk"])
def test_golden_gate(vector_fname):
    flank_left, flank_right = _golden_gate(vector_fname, "Esp3I", overhangs=("CACC", "GTTT"))
    assert len(flank_right) - len(flank_left) <= 1, "cutsite should be in the middle of the backbone"
    assert flank_left.endswith("CACC")
    assert flank_right.startswith("GTTT")

TEST_CASES_VARS = "query_fname,vector_fname,reference_plasmid,reference_spacer,expected_success,kwargs"
TEST_CASES = [
    ("pSPS217.28-SPS106.7.ab1", "tests/test_data/sanger/pLentiCRISPRv2.gbk", "pSPS217.28", "GGGGCCACTAGGGACAGGAT", False, dict()), # AAVS1
    ("SPS120.27-SPS106.5.ab1", "tests/test_data/sanger/CROPseq-Guide-Puro.gbk", "SPS120.27", "GTCACCAATCCTGTCCCTAG", True, dict(flank_width=4)), # AAVS1
]

@pytest.mark.parametrize(TEST_CASES_VARS, TEST_CASES, ids=[x[0] for x in TEST_CASES])
def test_validate_cloning(query_fname, vector_fname, reference_plasmid, reference_spacer, expected_success, kwargs):
    df = validate_cloning(
        query_dir="tests/test_data/sanger",
        spacers=pd.Series({reference_plasmid: reference_spacer}),
        vector=vector_fname,
        enzyme="Esp3I",
        overhangs=("CACC", "GTTT"),
        **kwargs,
        )
    assert df.shape == (1, 3)
    assert df.loc[query_fname, "success"] == expected_success
    assert df.loc[query_fname, "reference plasmid"] == reference_plasmid
    assert df.loc[query_fname, "spacer"] == reference_spacer

def test_cli_help():
    subprocess.run(["python", "-m", "be_scan", "validate_cloning", "--help"], check=True)

@pytest.mark.parametrize(TEST_CASES_VARS, TEST_CASES, ids=[x[0] for x in TEST_CASES])
def test_cli(query_fname, vector_fname, reference_plasmid, reference_spacer, expected_success, kwargs):
    in_fname = "tests/test_data/sanger/spacers.csv"
    out_fname = "tests/test_data/sanger/spacers_results.csv"
    try:
        # Create input reference file in desired format
        in_df = pd.Series({reference_plasmid: reference_spacer})
        in_df.rename_axis(index="plasmid", inplace=True)
        in_df.rename("spacer", inplace=True)
        in_df.to_csv(in_fname)
        subprocess.run(["python", "-m", "be_scan", "validate_cloning", 
                        "tests/test_data/sanger", # query_dir
                        in_fname,
                        vector_fname,
                        out_fname,
                        "--enzyme", "Esp3I",
                        "--overhangs", "CACC", "GTTT",
                        *[f"--{k}={v}" for k, v in kwargs.items()]
                        ], check=True)
        assert os.path.exists(out_fname)
        out_df = pd.read_csv(out_fname, index_col="sanger trace fname")
        assert out_df.shape == (1, 3)
        assert out_df.loc[query_fname, "success"] == expected_success
        assert out_df.loc[query_fname, "reference plasmid"] == reference_plasmid
        assert out_df.loc[query_fname, "spacer"] == reference_spacer
    finally:
        if os.path.exists(in_fname):
            os.remove(in_fname)
        if os.path.exists(out_fname):
            os.remove(out_fname)
