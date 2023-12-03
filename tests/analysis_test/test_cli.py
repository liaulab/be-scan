import subprocess

def test_entrypoint():
    out = subprocess.run("python -m be_scan -h", shell=True, capture_output=True)
    assert b"usage: be_scan" in out.stdout
    assert b"count_reads" in out.stdout
    assert b"average_reps" in out.stdout
    assert b"merge_and_norm" in out.stdout

# analysis CLI

def test_count_reads():
    out = subprocess.run("python -m be_scan count_reads -h", shell=True, capture_output=True)
    assert b"usage: be_scan count_reads" in out.stdout
    assert b"count the reads in the" in out.stdout

def test_average_reps():
    out = subprocess.run("python -m be_scan average_reps -h", shell=True, capture_output=True)
    assert b"usage: be_scan average_reps" in out.stdout
    assert b"Averages the replicates for each condition" in out.stdout

def test_merge_and_norm():
    out = subprocess.run("python -m be_scan merge_and_norm -h", shell=True, capture_output=True)
    assert b"usage: be_scan merge_and_norm" in out.stdout
    assert b"aggregate them into a single dataframe" in out.stdout
