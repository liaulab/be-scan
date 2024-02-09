import pytest

import pandas as pd
from pathlib import Path

from be_scan.sgrna.dataframes import *

testdata_filepath = "tests/test_data/sgrna/"
df1 = Path(f"{testdata_filepath}simpledf1.csv")
df2 = Path(f"{testdata_filepath}simpledf2.csv")
df3 = Path(f"{testdata_filepath}simpledf3.csv")

def test_merge_guide_df_pos(): 
    out = merge_guide_df(
        guide_df1_filepath =df1, 
        guide_df2_filepath =df2, 
        shared_col_names   =['sgRNA_seq', 'coding_seq'],
        sort_by            =['sgRNA_seq'],
        return_df          =True, 
        save_df            =False, 
    )
    assert out.shape[0] == 5 and out.shape[1] == 4

def test_add_guide_df_pos(): 
    out = add_guide_df(
        guides_df_filepath     =df1, 
        additional_df_filepath =df3, 
        return_df              =True, 
        save_df                =False, 
    )
    assert out.shape[0] == 4 and out.shape[1] == 3

# merging dataframes that don't have right sort_by column names inputs
def test_merge_guide_df_neg(): 
    with pytest.raises(AssertionError): 
        out = merge_guide_df(
            guide_df1_filepath =df1, 
            guide_df2_filepath =df2, 
            shared_col_names   =['sgRNA_seq', 'coding_seq'],
            sort_by            =['seq'],
            return_df          =True, 
            save_df            =False, 
        )

# adding dataframes without matching column names
def test_add_guide_df_neg(): 
    with pytest.raises(AssertionError): 
        out = add_guide_df(
            guides_df_filepath     =df1, 
            additional_df_filepath =df2, 
            return_df              =True, 
            save_df                =False, 
        )
