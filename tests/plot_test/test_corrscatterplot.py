import be_scan.plot as bsplot
import os
import pytest

file = "tests/test_data/plot_data/NZL10196_v9_comparisons.csv"

def check_comparisons(): 
    path = f"d3-posd3-neg_correlation_scatterplot.png"
    assert os.path.exists(path)
    os.remove(path)

# positive tests

def test_plot_corr_scatterplot_basic_pos():
    bsplot.plot_corr_scatterplot(df_filepath = file, condition1='d3-pos', condition2='d3-neg', 
    )
    check_comparisons()

def test_plot_corr_scatterplot_hue_pos():
    bsplot.plot_corr_scatterplot(df_filepath = file, condition1='d3-pos', condition2='d3-neg', 
                                 include_hue=True, hue_col='Mut_type', 
    )
    check_comparisons()
        
def test_plot_corr_scatterplot_filtervals_pos():
    bsplot.plot_corr_scatterplot(df_filepath = file, condition1='d3-pos', condition2='d3-neg', 
                                 filter_val=True, val_min=0.0, val_cols=['d3-pos', 'd6-pos', 'd9-pos'],
    )
    check_comparisons()

def test_plot_corr_scatterplot_filterparams_pos():
    bsplot.plot_corr_scatterplot(df_filepath = file, condition1='d3-pos', condition2='d3-neg', 
                                 filter_params=True, params_cols=['sgRNA_strand'], params_conditions=[['sense']],
    )
    check_comparisons()

# negative tests

@pytest.mark.parametrize("params", [{'hue_col':3}, # hue_col type error
                                    {'hue_col':'a'}, # hue_col not a column
                        ])
def test_plot_corr_scatterplot_hue_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.plot_corr_scatterplot(df_filepath = file, condition1='d3-pos', condition2='d3-neg', 
                                     include_hue=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{'val_min':'string', 'val_cols':['d3-pos', 'd6-pos', 'd9-pos']}, # val_min type error
                                    {'val_min':0.0, 'val_cols':'d3-pos'}, # val_cols type error
                        ])
def test_plot_corr_scatterplot_filtervals_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.plot_corr_scatterplot(df_filepath = file, condition1='d3-pos', condition2='d3-neg', 
                                     filter_val=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{'params_cols':'sgRNA_strand', 'params_conditions':[['sense']]}, # params_cols type error
                                    {'params_cols':['sgRNA_strand'], 'params_conditions':'sense'}, # params_conditions type error
                        ])
def test_plot_corr_scatterplot_filterparams_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.plot_corr_scatterplot(df_filepath = file, condition1='d3-pos', condition2='d3-neg', 
                                     filter_params=True, **params )
        check_comparisons()
