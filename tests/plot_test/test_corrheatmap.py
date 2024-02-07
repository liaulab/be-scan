import be_scan.plot
import os

file_dir = "tests/test_data/plot_data/"

def test_plot_corr_heatmap(monkeypatch):
    comparisons = ['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg']
    be_scan.plot.plot_corr_heatmap(df_filepath    =file_dir+'NZL10196_v9_comparisons.csv', 
                                   comparisons       =comparisons, 
    )
    path = "correlation_heatmap.pdf"
    assert os.path.exists(path)
    os.remove(path)
    