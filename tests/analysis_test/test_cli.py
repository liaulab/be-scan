import subprocess

def test_entrypoint():
    out = subprocess.run("python -m be_scan -h", shell=True, capture_output=True)
    assert b"usage: be_scan" in out.stdout
    assert b"count_reads" in out.stdout

def test_count_reads():
    out = subprocess.run("python -m be_scan count_reads -h", shell=True, capture_output=True)
    assert b"usage: be_scan count_reads" in out.stdout
    assert b"Count the reads in a FASTQ file" in out.stdout
