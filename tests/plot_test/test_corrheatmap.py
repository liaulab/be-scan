import be_scan.plot as bsplot
import os
import pytest

file = "tests/test_data/plot/NZL10196_v9_comparisons.csv"
comparisons = ['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg']
corr_params = {"df_filepath":file, "comparisons":comparisons, }

def check_comparisons(): 
    path = f"correlation_heatmap.png"
    assert os.path.exists(path)
    os.remove(path)

# positive tests

def test_plot_corrheatmap_basic_pos():
    bsplot.corr_heatmap(**corr_params, )
    check_comparisons()
        
def test_plot_corrheatmap_filtervals_pos():
    bsplot.corr_heatmap(**corr_params, filter_val=True, 
                             val_min=0.0, val_cols=['d3-pos', 'd6-pos', 'd9-pos'], )
    check_comparisons()

def test_plot_corrheatmap_filterparams_pos():
    bsplot.corr_heatmap(**corr_params, filter_params=True, 
                             params_cols=['sgRNA_strand'], params_conditions=[['sense']], )
    check_comparisons()

# negative tests

@pytest.mark.parametrize("params", [{'val_min':'string', 'val_cols':['d3-pos', 'd6-pos', 'd9-pos']}, # val_min type error
                                    {'val_min':0.0, 'val_cols':'d3-pos'}, # val_cols type error
                        ])
def test_plot_corrheatmap_filtervals_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.corr_heatmap(**corr_params, filter_val=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{'params_cols':'sgRNA_strand', 'params_conditions':[['sense']]}, # params_cols type error
                                    {'params_cols':['sgRNA_strand'], 'params_conditions':'sense'}, # params_conditions type error
                        ])
def test_plot_corrheatmap_filterparams_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.corr_heatmap(**corr_params, filter_params=True, **params )
        check_comparisons()
