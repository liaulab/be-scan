import subprocess

def test_entrypoint():
    out = subprocess.run("python -m be_scan -h", shell=True, capture_output=True)
    assert b"usage: be_scan" in out.stdout
    assert b"plot_scatterplot" in out.stdout
    assert b"plot_boxes" in out.stdout
    assert b"plot_corr_heatmap" in out.stdout
    assert b"plot_corr_scatterplot" in out.stdout

# plot CLI

def test_plot_scatterplot():
    out = subprocess.run("python -m be_scan plot_scatterplot -h", shell=True, capture_output=True)
    assert b"usage: be_scan plot_scatterplot" in out.stdout
    assert b"plots the data for each condition" in out.stdout

def test_plot_boxes():
    out = subprocess.run("python -m be_scan plot_boxes -h", shell=True, capture_output=True)
    assert b"usage: be_scan plot_boxes" in out.stdout
    assert b"plots chosen (ex control) guides by plot_column categories" in out.stdout

def test_plot_corr_heatmap():
    out = subprocess.run("python -m be_scan plot_corr_heatmap -h", shell=True, capture_output=True)
    assert b"usage: be_scan plot_corr_heatmap" in out.stdout
    assert b"scatterplot showing correlation between two given conditions" in out.stdout

def test_plot_corr_scatterplot():
    out = subprocess.run("python -m be_scan plot_corr_scatterplot -h", shell=True, capture_output=True)
    assert b"usage: be_scan plot_corr_scatterplot" in out.stdout
    assert b"heatmap showing correlation between all given comparison conditions" in out.stdout
