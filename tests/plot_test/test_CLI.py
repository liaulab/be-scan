import subprocess as sp
import pytest
import os

def test_entrypoint():
    out = sp.run("python -m be_scan -h", shell=True, capture_output=True)
    assert b"usage: be_scan" in out.stdout
    assert b"plot_scatterplot" in out.stdout
    assert b"boxplot" in out.stdout
    assert b"corr_heatmap" in out.stdout
    assert b"corr_jointplot" in out.stdout

# plot CLI

def test_plot_scatterplot():
    out = sp.run("python -m be_scan plot_scatterplot -h", shell=True, capture_output=True)
    assert b"usage: be_scan plot_scatterplot" in out.stdout
    assert b"plots the data for each condition" in out.stdout

def test_boxplot():
    out = sp.run("python -m be_scan boxplot -h", shell=True, capture_output=True)
    assert b"usage: be_scan boxplot" in out.stdout
    assert b"plots chosen (ex control) guides by plot_column categories" in out.stdout

def test_corr_heatmap():
    out = sp.run("python -m be_scan corr_heatmap -h", shell=True, capture_output=True)
    assert b"usage: be_scan corr_heatmap" in out.stdout
    assert b"scatterplot showing correlation between two given conditions" in out.stdout

def test_corr_jointplot():
    out = sp.run("python -m be_scan corr_jointplot -h", shell=True, capture_output=True)
    assert b"usage: be_scan corr_jointplot" in out.stdout
    assert b"heatmap showing correlation between all given comparison conditions" in out.stdout

# more specific tests
    
plot_file = "tests/test_data/plot/NZL10196_v9_comparisons.csv"

boxes_query1 = "Domain -c d3-pos d3-neg d6-pos -pc PWWP ADD MTase Nterm" 
boxes_query2 = "--savefig --show --out_name boxes.png"
@pytest.mark.parametrize("query", [
    "", 
    "--neg_ctrl --neg_ctrl_col Gene --neg_ctrl_conditions NON-GENE", 
    "--filter_val --val_min 0.0 --val_cols d3-pos d6-pos d9-pos", 
    ])
def test_boxplot_integration_pos(query): 
    out = sp.run(f"python3 -m be_scan boxplot {plot_file} {boxes_query1} {query} {boxes_query2}",
                 shell=True, capture_output=True)
    assert out.returncode == 0


scatter_query1 = "Edit_site_3A1 -c d3-pos d3-neg d6-pos" 
scatter_query2 = "--savefig --show --out_name scatterplot.png"
@pytest.mark.parametrize("query", [
    "", 
    "--include_hue --hue_col Mut_type", 
    "--neg_ctrl --neg_ctrl_col Gene --neg_ctrl_conditions NON-GENE", 
    "--filter_val --val_min 0.0 --val_cols d3-pos d6-pos d9-pos", 
    "--autoannot --autoannot_label sgRNA_ID --autoannot_top 10", 
    "--autoannot --autoannot_label sgRNA_ID --autoannot_cutoff 2.0", 
    ])
def test_plot_scatterplot_integration_pos(query): 
    out = sp.run(f"python3 -m be_scan plot_scatterplot {plot_file} {scatter_query1} {query} {scatter_query2}",
                 shell=True, capture_output=True)
    assert out.returncode == 0


corrscatter_query1 = "d3-pos d6-pos" 
corrscatter_query2 = "--savefig --show --out_name corrscatterplot.png"
@pytest.mark.parametrize("query", [
    "", 
    "--include_hue --hue_col Mut_type", 
    "--filter_val --val_min 0.0 --val_cols d3-pos d6-pos d9-pos", 
    ])
def test_plot_corrscatterplot_integration_pos(query): 
    out = sp.run(f"python3 -m be_scan corr_jointplot {plot_file} {corrscatter_query1} {query} {corrscatter_query2}",
                 shell=True, capture_output=True)
    assert out.returncode == 0


corrheat_query1 = "-c d3-pos d6-pos d9-pos"
corrheat_query2 = "--savefig --show --out_name corrheatmap.png"
@pytest.mark.parametrize("query", [
    "", 
    "--filter_val --val_min 0.0 --val_cols d3-pos d6-pos d9-pos", 
    ])
def test_plot_corrheatmap_integration_pos(query): 
    out = sp.run(f"python3 -m be_scan corr_heatmap {plot_file} {corrheat_query1} {query} {corrheat_query2}",
                 shell=True, capture_output=True)
    assert out.returncode == 0
