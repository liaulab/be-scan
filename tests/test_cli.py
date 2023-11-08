import subprocess

def test_entrypoint():
    out = subprocess.run("python -m bigscam -h", shell=True, capture_output=True)
    assert b"usage: bigscam" in out.stdout
    assert b"count_reads" in out.stdout

def test_count_reads():
    out = subprocess.run("python -m bigscam count_reads -h", shell=True, capture_output=True)
    assert b"usage: bigscam count_reads" in out.stdout
    assert b"Count the reads in a FASTQ file" in out.stdout
