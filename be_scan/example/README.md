## BE-scan example data

be_scan/example/SampleData-PRC2-[WT-Dropout-ABE].csv contains a data table with an annotated sgRNA library and sgRNA counts for three replicates of multiple treatments and timepoints. These data can be used to demo BE-scan's analysis features. To use, read the data into a pandas dataframe. Then run, in order:
1. be_scan.analysis.log_transform() to merge replicates of different samples and log2 transform data, and output this as a csv file.
2. be_scan.analysis.compare_conds() to calculate log fold change between conditions, and output this as a csv file.
3. be_scan.calc_controls() to z-score each each sgRNA relative to the negative control sgRNAs for that comparison.
Each of these functions take less than 2 minutes to run.
More details on each function are available at https://be-scan.readthedocs.io/en/main/index.html. 
