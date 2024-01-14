import be_scan.plot
import os

file_dir = "tests/test_data/plot_data/"

def test_plot_scatterplot(monkeypatch):
    comparisons = ['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg']
    be_scan.plot.plot_scatterplot(df_filepath       =file_dir+'NZL10196_v9_comparisons.csv', 
                                  x_column          ='Edit_site_3A1', 
                                  comparisons       =comparisons, 
                                  hue=True, hue_column='Mut_type', 
                                  neg_ctrl=True, neg_ctrl_col='Gene', neg_ctrl_conditions=['NON-GENE'],
    )
    for suffix in comparisons: 
        path = f"scatterplot{suffix}.pdf"
        assert os.path.exists(path)
        os.remove(path)

def test_plot_corr_heatmap(monkeypatch):
    comparisons = ['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg']
    be_scan.plot.plot_corr_heatmap(df_filepath    =file_dir+'NZL10196_v9_comparisons.csv', 
                                   comparisons       =comparisons, 
    )
    path = "correlation_heatmap.pdf"
    assert os.path.exists(path)
    os.remove(path)

def test_plot_corr_scatterplot(monkeypatch):
    c1 = 'd3-neg'
    c2 = 'd9-pos'
    be_scan.plot.plot_corr_scatterplot(df_filepath =file_dir+'NZL10196_v9_comparisons.csv', 
                                       condition1  =c1, 
                                       condition2  =c2, 
                                       hue_column  ='Mut_type',
    )
    path = f"{c1}{c2}_correlation_scatterplot.pdf"
    assert os.path.exists(path)
    os.remove(path)

def test_plot_boxes(monkeypatch):
    comparisons = ['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg']
    be_scan.plot.plot_boxes(df_filepath       =file_dir+'NZL10196_v9_comparisons.csv', 
                            plot_column       ='Domain', 
                            plot_conditions   =['PWWP', 'ADD', 'MTase', 'Nterm'], 
                            comparisons       =comparisons, 
                            neg_ctrl=True, neg_ctrl_col='Gene', neg_ctrl_conditions=['NON-GENE'],
    )
    for suffix in comparisons: 
        path = f"boxes{suffix}.pdf"
        assert os.path.exists(path)
        os.remove(path)
