import be_scan.plot
import os

file_dir = "tests/test_data/plot_data/"

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
    