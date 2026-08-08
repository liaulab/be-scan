import be_scan.plot as bsplot
import os
import pytest
import pandas as pd

df_ezh2 = pd.read_csv("tests/test_data/plot/q575r_ezh2_abe_missense_pos.csv")
af_ezh2 = "tests/test_data/plot/AF_Q15910.pdb"
pdb_ezh2 = "tests/test_data/plot/6wkr-chainC-EZH2.pdb"

df_prc2 = pd.read_csv("tests/test_data/plot/q575r_abe_pos.csv")
af_prc2 = "tests/test_data/plot/AF-PRC2-EZH2iso2.pdb"
pdb_prc2 = "tests/test_data/plot/6wkr.pdb"

pwes_params = {
    "x_col": "xpos", "scores_col": "DMSO-PRC2i", 
    "gauss_std": 16, "dend_t": 14, "tanh_a": 1.0, 
    "out_prefix": "test", 
    "out_dir": "tests/test_data/plot/test", 
    }

files = [
    'test_cluster_boxplots.pdf', 'test_cluster_boxplots.png', 
    'test_cluster_heatmap.pdf', 'test_cluster_heatmap.png', 
    'test_cluster_histogram.pdf', 'test_cluster_histogram.png', 
    'test_combined_legend.pdf', 'test_combined_legend.png', 
    'test_PWES_heatmap.pdf', 'test_PWES_heatmap.png', 
    'test_scatterplot_subplots.pdf', 'test_scatterplot_subplots.png', 
    'test_scatterplot.pdf', 'test_scatterplot.png', 
    'test_scatterplotlegend.pdf', 'test_scatterplotlegend.png', 
    'test_aas_dict.pickle', 
]

def remove_outputs(): 
    for f in files: 
        path = f"tests/test_data/plot/test/{f}"
        if os.path.exists(path): 
            print('Deleting', path)
            os.remove(path)

# positive tests

# pdb_file = EZH2_pdb, 
# df_scores = pd.read_csv(files_whole[0]), 
# domains_list=EZH2_domains_list, 

def test_clustering_ezh2_af(): 
    temp = bsplot.pwes_clustering(
        **pwes_params, 
        pdb_file = af_ezh2, 
        df_scores = df_ezh2, 
        show_plots=False, 
        )
    del temp
    remove_outputs()

def test_clustering_ezh2_pdb(): 
    temp = bsplot.pwes_clustering(
        **pwes_params, 
        pdb_file = pdb_ezh2, 
        df_scores = df_ezh2, 
        show_plots=False, 
        )
    del temp
    remove_outputs()

def test_clustering_prc2_af(): 
    temp = bsplot.pwes_clustering(
        **pwes_params, 
        pdb_file = af_prc2, 
        df_scores = df_prc2, 
        show_plots=False, 
        )
    del temp
    remove_outputs()

def test_clustering_prc2_pdb(): 
    temp = bsplot.pwes_clustering(
        **pwes_params, 
        pdb_file = pdb_prc2, 
        df_scores = df_prc2, 
        show_plots=False, 
        )
    del temp
    remove_outputs()

# def test_clustering_prc2_af_domains(): 
#     temp = bsplot.pwes_clustering(
#         **pwes_params, 
#         pdb_file = af_prc2, 
#         df_scores = df_prc2, 
#         show_plots=False, 
#         )
#     del temp
#     remove_outputs()

# def test_clustering_prc2_pdb_domains(): 
#     temp = bsplot.pwes_clustering(
#         **pwes_params, 
#         pdb_file = pdb_prc2, 
#         df_scores = df_prc2, 
#         show_plots=False, 
#         )
#     del temp
#     remove_outputs()

# negative tests

# @pytest.mark.parametrize("params", [{"hue_col":3}, # hue_col type error
#                                     {"hue_col":"a"}, # hue_col not a column
#                         ])
# def test_scatterplot_hue_neg(params):
#     with pytest.raises(AssertionError): 
#         bsplot.scatterplot(**scatter_params,
#                                 include_hue=True, **params )
#         check_comparisons()

# @pytest.mark.parametrize("params", [{"neg_ctrl_col":["Gene"], "neg_ctrl_conditions":["NON-GENE"]}, # neg_ctrl_col type error
#                                     {"neg_ctrl_col":"Gene", "neg_ctrl_conditions":"NON-GENE"}, # neg_ctrl_conditions type error
#                                     {"neg_ctrl_col":"a", "neg_ctrl_conditions":["NON-GENE"]}, # neg_ctrl_col not a column
#                         ])
# def test_scatterplot_negctrl_neg(params):
#     with pytest.raises(AssertionError): 
#         bsplot.scatterplot(**scatter_params, neg_ctrl=True, **params )
#         check_comparisons()

# @pytest.mark.parametrize("params", [{"val_min":"string", "val_cols":["d3-pos", "d6-pos", "d9-pos"]}, # val_min type error
#                                     {"val_min":0, "val_cols":"d3-pos"}, # val_cols type error
#                         ])
# def test_scatterplot_filtervals_neg(params):
#     with pytest.raises(AssertionError): 
#         bsplot.scatterplot(**scatter_params, filter_val=True, **params )
#         check_comparisons()

# @pytest.mark.parametrize("params", [{"params_cols":"sgRNA_strand", "params_conditions":[["sense"]]}, # params_cols type error
#                                     {"params_cols":["sgRNA_strand"], "params_conditions":"sense"}, # params_conditions type error
#                         ])
# def test_scatterplot_filterparams_neg(params):
#     with pytest.raises(AssertionError): 
#         bsplot.scatterplot(**scatter_params, filter_params=True, **params )
#         check_comparisons()

# @pytest.mark.parametrize("params", [{"autoannot_label":10, "autoannot_top":10}, # autoannot_label type error
#                                     {"autoannot_label":"a", "autoannot_top":10}, # autoannot_label not a column error
#                                     {"autoannot_label":"sgRNA_ID", "autoannot_top":"a"}, # autoannot_top type error
#                         ])
# def test_scatterplot_autoannot1_neg(params):
#     with pytest.raises(AssertionError): 
#         bsplot.scatterplot(**scatter_params, autoannot=True, **params )
#         check_comparisons()

# @pytest.mark.parametrize("params", [{"autoannot_label":10, "autoannot_cutoff":1.0}, # autoannot_label type error
#                                     {"autoannot_label":"a", "autoannot_cutoff":1.0}, # autoannot_label not a column error
#                                     {"autoannot_label":"sgRNA_ID", "autoannot_cutoff":'x'}, # autoannot_cutoff type error
#                         ])
# def test_scatterplot_autoannot2_neg(params):
#     with pytest.raises(AssertionError): 
#         bsplot.scatterplot(**scatter_params, autoannot=True, **params )
#         check_comparisons()
