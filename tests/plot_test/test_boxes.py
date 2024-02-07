import be_scan.plot
import os

file_dir = "tests/test_data/plot_data/"


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
