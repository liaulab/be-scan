import subprocess
import os

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


file_dir = "tests/test_data/analysis/"


sample_sheet1 = f"{file_dir}guides_sample_sheet.csv"
annotated_lib1 = f"{file_dir}guides_ref.csv"
def test_count_reads_integration_pos(): 
    out = subprocess.run(f"python3 -m be_scan count_reads {sample_sheet1} {annotated_lib1} --file_dir {file_dir} --out_dir {file_dir}", 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm {0}".format(f"{file_dir}counts_library.csv"))


sample_sheet2 = f"{file_dir}sample_sheet.csv"
annotated_lib2 = f"{file_dir}count_reads_sample_out.csv"
def test_merge_and_norm_integration_pos(): 
    out = subprocess.run(f"python3 -m be_scan merge_and_norm {sample_sheet2} {annotated_lib2} --out_dir {file_dir} --controls counts1", 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm {0}".format(f"{file_dir}agg_log2_t0.csv"))


annotated_lib3 = f"{file_dir}merge_and_norm_sample_out.csv"
def test_average_reps_integration_pos(): 
    out = subprocess.run(f"python3 -m be_scan average_reps {sample_sheet2} {annotated_lib3} --out_dir {file_dir}", 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm {0}".format(f"{file_dir}avg_conds.csv"))


comparisons = f"{file_dir}comparisons.csv"
annotated_lib4 = f"{file_dir}compare_conds_sample_out.csv"
def test_compare_conds_integration_pos(): 
    out = subprocess.run(f"python3 -m be_scan compare_conds {comparisons} {annotated_lib4} --out_dir {file_dir}", 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm {0}".format(f"{file_dir}conditions.csv"))


annotated_lib5 = f"{file_dir}compare_conds_sample_out.csv"
def test_calc_controls_integration_pos(): 
    out = subprocess.run(f"python3 -m be_scan calc_controls {annotated_lib5} gene -sc cond1 -ncc control --out_dir {file_dir}", 
                         shell=True, capture_output=True)
    assert out.returncode == 0
    os.system("rm {0}".format(f"{file_dir}stats.txt"))