import be_scan.plot as bsplot
import os
import pytest

file = "tests/test_data/plot_data/NZL10196_v9_comparisons.csv"
comparisons = ['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg']
plot_conditions = ['PWWP', 'ADD', 'MTase', 'Nterm']

def check_comparisons(): 
    for suffix in comparisons: 
        path = f"boxes{suffix}.png"
        assert os.path.exists(path)
        os.remove(path)

# positive tests
        
def test_plot_boxes_basic_pos():
    bsplot.plot_boxes(df_filepath = file, plot_column = 'Domain', 
                      comparisons = comparisons, plot_conditions = plot_conditions,
    )
    check_comparisons()

def test_plot_boxes_negctrl_pos():
    bsplot.plot_boxes(df_filepath = file, plot_column = 'Domain', 
                      comparisons = comparisons, plot_conditions = plot_conditions,
                      neg_ctrl=True, neg_ctrl_col='Gene', neg_ctrl_conditions=['NON-GENE'],
    )
    check_comparisons()
        
def test_plot_boxes_filtervals_pos():
    bsplot.plot_boxes(df_filepath = file, plot_column = 'Domain', 
                      comparisons = comparisons, plot_conditions = plot_conditions,
                      filter_val=True, val_min=0.0, val_cols=['d3-pos', 'd6-pos', 'd9-pos'],
    )
    check_comparisons()

def test_plot_boxes_filterparams_pos():
    bsplot.plot_boxes(df_filepath = file, plot_column = 'Domain', 
                      comparisons = comparisons, plot_conditions = plot_conditions,
                      filter_params=True, params_cols=['sgRNA_strand'], params_conditions=[['sense']],
    )
    check_comparisons()

# negative tests

@pytest.mark.parametrize("params", [{'neg_ctrl_col':['Gene'], 'neg_ctrl_conditions':['NON-GENE']}, # neg_ctrl_col type error
                                    {'neg_ctrl_col':'Gene', 'neg_ctrl_conditions':'NON-GENE'}, # neg_ctrl_conditions type error
                                    {'neg_ctrl_col':'a', 'neg_ctrl_conditions':['NON-GENE']}, # neg_ctrl_col not a column
                        ])
def test_plot_boxes_negctrl_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.plot_boxes(df_filepath = file, plot_column = 'Domain', 
                          comparisons = comparisons, plot_conditions = plot_conditions,
                          neg_ctrl=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{'val_min':'string', 'val_cols':['d3-pos', 'd6-pos', 'd9-pos']}, # val_min type error
                                    {'val_min':0, 'val_cols':'d3-pos'}, # val_cols type error
                        ])
def test_plot_boxes_filtervals_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.plot_boxes(df_filepath = file, plot_column = 'Domain', 
                          comparisons = comparisons, plot_conditions = plot_conditions,
                          filter_val=True, **params )
        check_comparisons()

@pytest.mark.parametrize("params", [{'params_cols':'sgRNA_strand', 'params_conditions':[['sense']]}, # params_cols type error
                                    {'params_cols':['sgRNA_strand'], 'params_conditions':'sense'}, # params_conditions type error
                        ])
def test_plot_boxes_filterparams_neg(params):
    with pytest.raises(AssertionError): 
        bsplot.plot_boxes(df_filepath = file, plot_column = 'Domain', 
                          comparisons = comparisons, plot_conditions = plot_conditions,
                          filter_params=True, **params )
        check_comparisons()

