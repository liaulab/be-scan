import be_scan.plot as bsplot
import os
import pytest

file = "tests/test_data/plot/NZL10196_v9_comparisons.csv"
comparisons = ["d3-pos", "d3-neg", "d6-pos", "d6-neg", "d9-pos", "d9-neg"]
scatter_params = {"df_filepath":file, "x_column":"Edit_site_3A1", "comparisons":comparisons}

def check_comparisons(): 
    for suffix in comparisons: 
        path = f"scatterplot{suffix}.png"
        assert os.path.exists(path)
        os.remove(path)

# positive tests
        
def test_scatterplot_basic_pos():
    bsplot.scatterplot(**scatter_params, )
    check_comparisons()

def test_scatterplot_hue_pos():
    bsplot.scatterplot(**scatter_params,
                            include_hue=True, hue_col="Mut_type", )
    check_comparisons()

def test_scatterplot_negctrl_pos():
    bsplot.scatterplot(**scatter_params, neg_ctrl=True, 
                            neg_ctrl_col="Gene", neg_ctrl_conditions=["NON-GENE"], )
    check_comparisons()
        
def test_scatterplot_filtervals_pos():
    bsplot.scatterplot(**scatter_params, filter_val=True, 
                            val_min=0.0, val_cols=["d3-pos", "d6-pos", "d9-pos"], )
    check_comparisons()

def test_scatterplot_filterparams_pos():
    bsplot.scatterplot(**scatter_params, filter_params=True, 
                            params_cols=["sgRNA_strand"], params_conditions=[["sense"]], )
    check_comparisons()

def test_scatterplot_autoannot1_pos():
    bsplot.scatterplot(**scatter_params, autoannot=True, 
                            autoannot_label="sgRNA_ID", autoannot_top=10, )
    check_comparisons()

def test_scatterplot_autoannot2_pos():
    bsplot.scatterplot(**scatter_params, autoannot=True, 
                            autoannot_label="sgRNA_ID", autoannot_cutoff=2.0, )
    check_comparisons()

# negative tests

@pytest.mark.parametrize("params", [{"hue_col":3}, # hue_col type error
                                    {"hue_col":"a"}, # hue_col not a column
                        ])
def test_scatterplot_hue_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.scatterplot(**scatter_params,
                                include_hue=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{"neg_ctrl_col":["Gene"], "neg_ctrl_conditions":["NON-GENE"]}, # neg_ctrl_col type error
                                    {"neg_ctrl_col":"Gene", "neg_ctrl_conditions":"NON-GENE"}, # neg_ctrl_conditions type error
                                    {"neg_ctrl_col":"a", "neg_ctrl_conditions":["NON-GENE"]}, # neg_ctrl_col not a column
                        ])
def test_scatterplot_negctrl_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.scatterplot(**scatter_params, neg_ctrl=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{"val_min":"string", "val_cols":["d3-pos", "d6-pos", "d9-pos"]}, # val_min type error
                                    {"val_min":0, "val_cols":"d3-pos"}, # val_cols type error
                        ])
def test_scatterplot_filtervals_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.scatterplot(**scatter_params, filter_val=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{"params_cols":"sgRNA_strand", "params_conditions":[["sense"]]}, # params_cols type error
                                    {"params_cols":["sgRNA_strand"], "params_conditions":"sense"}, # params_conditions type error
                        ])
def test_scatterplot_filterparams_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.scatterplot(**scatter_params, filter_params=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{"autoannot_label":10, "autoannot_top":10}, # autoannot_label type error
                                    {"autoannot_label":"a", "autoannot_top":10}, # autoannot_label not a column error
                                    {"autoannot_label":"sgRNA_ID", "autoannot_top":"a"}, # autoannot_top type error
                        ])
def test_scatterplot_autoannot1_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.scatterplot(**scatter_params, autoannot=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{"autoannot_label":10, "autoannot_cutoff":1.0}, # autoannot_label type error
                                    {"autoannot_label":"a", "autoannot_cutoff":1.0}, # autoannot_label not a column error
                                    {"autoannot_label":"sgRNA_ID", "autoannot_cutoff":'x'}, # autoannot_cutoff type error
                        ])
def test_scatterplot_autoannot2_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.scatterplot(**scatter_params, autoannot=True, **params )
        check_comparisons()
