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
        