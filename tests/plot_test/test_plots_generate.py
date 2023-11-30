import be_scan.plot

def test_plot_scatterplot(monkeypatch):
    be_scan.plot.plot_scatterplot(df_filepath     ='tests/test_data/plot_test_data/NZL10196_v9_comparisons.csv', 
                                x_column          ='Edit_site_3A1', 
                                y_column          ='log2_fc', 
                                hue_column        ='Mut_type', 
                                comparisons       =['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg'], 
                                neg_ctrl_col      ='Gene', 
                                neg_ctrl_category ='NON-GENE',
                                xmin              =200, 
                                savefig           =False,
                                )
    assert True

def test_plot_corr_heatma(monkeypatch):
    be_scan.plot.plot_corr_heatmap(df_filepath    ='tests/test_data/plot_test_data/NZL10196_v9_comparisons.csv', 
                                comparisons       =['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg'], 
                                
                                savefig           =False,
                                )
    assert True

def test_plot_corr_scatterplot(monkeypatch):
    be_scan.plot.plot_corr_scatterplot(df_filepath ='tests/test_data/plot_test_data/NZL10196_v9_comparisons.csv', 
                                   condition1  ='d3-neg', 
                                   condition2  ='d9-pos', 
                                   hue_column  ='Mut_type',
                                   
                                   savefig     =False,
                                   )
    assert True

def test_plot_boxes(monkeypatch):
    be_scan.plot.plot_boxes(df_filepath       ='tests/test_data/plot_test_data/NZL10196_v9_comparisons.csv', 
                            plot_column       ='Domain', 
                            plot_conditions   =['PWWP', 'ADD', 'MTase'], 
                            y_column          ='log2_fc', 
                            comparisons       =['d3-pos', 'd3-neg', 'd6-pos', 'd6-neg', 'd9-pos', 'd9-neg'], 
                            neg_ctrl_col      ='Gene', 
                            neg_ctrl_category ='NON-GENE',

                            savefig           =False,
    )
    assert True
