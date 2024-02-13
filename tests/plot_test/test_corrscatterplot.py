import be_scan.plot as bsplot
import os
import pytest

file = "tests/test_data/plot/NZL10196_v9_comparisons.csv"
corr_params = {"df_filepath":file, "condition1":"d3-pos", "condition2":"d3-neg", }

def check_comparisons(): 
    path = f"d3-posd3-neg_correlation_jointplot.png"
    assert os.path.exists(path)
    os.remove(path)

# positive tests

def test_corr_jointplot_basic_pos():
    bsplot.corr_jointplot(**corr_params, )
    check_comparisons()

def test_corr_jointplot_hue_pos():
    bsplot.corr_jointplot(**corr_params,
                                 include_hue=True, hue_col="Mut_type", )
    check_comparisons()
        
def test_corr_jointplot_filtervals_pos():
    bsplot.corr_jointplot(**corr_params, filter_val=True, 
                                 val_min=0.0, val_cols=["d3-pos", "d6-pos", "d9-pos"], )
    check_comparisons()

def test_corr_jointplot_filterparams_pos():
    bsplot.corr_jointplot(**corr_params, filter_params=True, 
                                 params_cols=["sgRNA_strand"], params_conditions=[["sense"]], )
    check_comparisons()

# negative tests

@pytest.mark.parametrize("params", [{"hue_col":3}, # hue_col type error
                                    {"hue_col":"a"}, # hue_col not a column
                        ])
def test_corr_jointplot_hue_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.corr_jointplot(**corr_params, include_hue=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{"val_min":"string", "val_cols":["d3-pos", "d6-pos", "d9-pos"]}, # val_min type error
                                    {"val_min":0.0, "val_cols":"d3-pos"}, # val_cols type error
                        ])
def test_corr_jointplot_filtervals_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.corr_jointplot(**corr_params, filter_val=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{"params_cols":"sgRNA_strand", "params_conditions":[["sense"]]}, # params_cols type error
                                    {"params_cols":["sgRNA_strand"], "params_conditions":"sense"}, # params_conditions type error
                        ])
def test_corr_jointplot_filterparams_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.corr_jointplot(**corr_params, filter_params=True, **params )
        check_comparisons()
