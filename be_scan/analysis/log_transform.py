"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Aggregate raw counts and perform log2-transform and t0 normalization}
"""

from pathlib import Path

import numpy as np
import pandas as pd

def log_transform(
    sample_sheet, 
    library_counts, 
                   
    controls=['t0'], 
    in_dir='', out_dir='', out_file='library_LFC.csv', 
    save=True, return_df=True,
    ):
    
    """[Summary]
    For a given set of samples and their raw read count files from count_reads,
    aggregate them into a single dataframe, normalize counts to reads per million,
    perform log2 transform, and then normalize to t0. The aggregated raw reads,
    log2 transformed values, and t0 normalized values can be saved to csv.

    Parameters
    ------------
    sample_sheet : str or path or pandas DataFrame()
        REQUIRED COLS: 'condition', 'counts_file'
        a sheet with information on sequence id, 
        in_fastq (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_np (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    library_counts : str or path or pandas DataFrame()
        String or path to the reference file. library_counts must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.

    controls : str, default ['t0']
        Name of the control condition samples in sample_sheet. 
    in_dir : str or path, defaults to ''
        String or path to the directory where library_counts and sample_sheet is found. 
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'agg_log2_t0.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe

    Returns
    ------------
    """
    in_path = Path(in_dir)
    out_path = Path(out_dir)

    # IMPORT FILES, CHECK FOR REQUIREMENTS #
    if isinstance(library_counts, (str, Path)):
        df_ref = pd.read_csv(in_path / library_counts)
    else: df_ref = library_counts.copy()
    if isinstance(sample_sheet, (str, Path)):
        df_samples = pd.read_csv(in_path / sample_sheet)
    else: df_samples = sample_sheet.copy()

    dict_counts = dict(zip(df_samples.condition, df_samples.counts_file))
    dict_conds = dict(zip(df_samples.condition, df_samples.agg_conditions))
    df_map = pd.DataFrame(data=dict_conds.items(), columns=['rep','condition'])

    # CHECK CONTROLS IN DF AND AVERAGE IF NECESSARY #
    conditions = df_samples.condition.tolist()
    for control in controls: 
        if control not in dict_counts.keys(): 
            raise Exception (f"{sample_sheet} is missing the {control} sample")
    if len(controls) > 0: 
        df_ref['control_avg'] = df_ref[controls].mean(axis=1)
        conditions = ['control_avg'] + conditions
    
    # CALCULATE LFC AND CONTROLS #
    for sample in conditions: 
        # log2 normalization
        total_reads = pd.to_numeric(df_ref[sample]).sum()
        suffix = '_LFC'
        df_ref[sample+suffix] = df_ref[sample].apply(lambda x: np.log2((x * 1000000 / total_reads) + 1))
        if sample not in controls and len(controls) > 0: 
            # t0 normalization
            suffix = '_LFCminusControl'
            df_ref[sample+suffix] = df_ref[sample+'_LFC'].sub(df_ref['control_avg_LFC'])

    # AVERAGE LFC #
    for cond in df_map['condition'].unique().tolist(): 
        if cond not in controls: 
            # for each unique condition, find the reps
            reps = [x+suffix for x in df_map.loc[df_map['condition'] == cond]['rep'].tolist()]
            # skip averaging for single replicates (otherwise breaks script)
            if len(reps) > 1: 
                df_ref[cond+suffix+'_avg'] = df_ref[reps].mean(axis=1)
                df_ref[cond+suffix+'_stdev'] = df_ref[reps].std(axis=1)
            elif len(reps) == 1: 
                df_ref[cond+suffix+'_avg'] = df_ref[reps]
                df_ref[cond+suffix+'_stdev'] = 0
            else: 
                raise Exception('Error! Replicate number not valid')

    # SAVE DF AND RETURN #
    if save == True: 
        Path.mkdir(out_path, exist_ok=True)
        df_ref.to_csv(out_path / out_file, index=False)
        print('merge_and_norm output to', str(out_path / out_file))
    print('Merge and normalize completed')
    if return_df: 
        return df_ref
    