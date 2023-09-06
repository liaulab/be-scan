# -*- coding: utf-8 -*-
"""
@author: Kevin Ngan

Code for analysis of CRISPR tiling screens, including functions for QC and
batch processing. See function docstrings for exhaustive documentation.
"""
#%% import packages

from collections import Counter
from pathlib import Path
import time
import warnings

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from Bio import SeqIO

#%% count_reads() - count sgRNA reads in FASTQ (adapted from count_spacers)

def count_reads(in_fastq, in_ref, KEY_INTERVAL=(10,80), DIR='FWD',
                KEY='CGAAACACCG', KEY_REV='GTTTTAGA', out_counts='counts.csv',
                out_np='np_counts.csv', out_stats='stats.txt'):
    """
    Count the reads in a FASTQ file and assign them to a reference sgRNA set.

    Given a set of sgRNA sequences and a FASTQ file, count the reads in the
    FASTQ, assign the reads to sgRNAs, and export the counts to a csv file. All
    sgRNA sequences not found in the reference file (non-perfect matches) are
    written to a separate csv file (npcounts).

    Parameters
    ----------
    in_fastq : str or path
        String or path to the FASTQ file to be processed.
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY : str, default 'CGAAACACCG'
        Upstream sequence that identifies the position of the sgRNA when used
        with DIR='FWD'. The default is the end of the hU6 promoter.
    KEY_REV : str, default 'GTTTTAGA'
        Downstream sequence that identifies the position of the sgRNA when used
        with DIR='REV'. The default is the start of the sgRNA scaffold sequence.
    DIR : {'FWD', 'REV'}, default 'FWD'
        The direction to identify the position of the sgRNA relative to the
        KEY sequence. 'FWD' is upstream of sgRNA, 'REV' is downstream of sgRNA.
    out_counts : str or path, default 'counts.csv'
        String or path for the output csv file with perfect sgRNA matches.
    out_np : str or path, default 'np_counts.csv'
        String or path for the output csv file with non-perfect sgRNA matches.
    out_stats : str or path, default 'stats.txt'
        String or path for the output txt file with the read counting statistics.
    """

    # STEP 1A: OPEN INPUT FILES FOR PROCESSING, CHECK FOR REQUIRED FORMATTING
    # look for 'sgRNA_seq' column, raise Exception if missing
    df_ref = pd.read_csv(in_ref, header=0) # explicit header = first row
    if 'sgRNA_seq' not in df_ref.columns.tolist():
        raise Exception('in_ref is missing column: sgRNA_seq')
    # look for other cols, raise Warning if suggested cols are missing
    list_headcols = ['sgRNA_ID', 'sgRNA_seq', 'Gene', 'cut_site_AA', 'Domain']
    if not all(col in df_ref.columns.tolist() for col in list_headcols):
        list_miss = [col for col in list_headcols if col not in df_ref.columns.tolist()]
        warnings.warn('Warning! in_ref is missing column(s) for downstream functions: ' + str(list_miss))
    # try opening input FASTQ, raise Exception if not possible
    try:
        handle = open(in_fastq)
    except:
        print('Error! Could not open the FASTQ file: %s' % in_fastq)
        return

    # STEP 1B: SET UP VARIABLES FOR SCRIPT
    # make dictionary to hold sgRNA counts - sgRNA_seq, count as k,v
    dict_perfects = {sgRNA:0 for sgRNA in df_ref['sgRNA_seq']}
    list_np = [] # placeholder list for non-perfect matches
    num_reads = 0 # total number of reads processed
    num_perfect_matches = 0 # count of reads with a perfect match to library
    num_np_matches = 0 # count of reads without a perfect match to library
    num_nokey = 0 # count of reads where key was not found
    KEY_START, KEY_END = KEY_INTERVAL[0], KEY_INTERVAL[1] # set the key interval

    # STEP 2: PROCESS FASTQ FILE READS AND ADD COUNTS TO DICT
    readiter = SeqIO.parse(handle, 'fastq') # process reads in fastq file
    # find sgRNA using FORWARD direction (default)
    if DIR == 'FWD':
        for record in readiter: # contains the seq and Qscore etc.
            num_reads += 1
            read_sequence = str.upper(str(record.seq))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY)
            if key_index >= 0: # if key found
                start_index = key_index + KEY_START + len(KEY)
                guide = read_sequence[start_index:(start_index + 20)]
                if guide in dict_perfects:
                    dict_perfects[guide] += 1
                    num_perfect_matches += 1
                else:
                    num_np_matches += 1
                    list_np.append(guide)
            else:
                num_nokey += 1
    # find sgRNA using REVERSE direction
    elif DIR == 'REV':
        for record in readiter: # contains the seq and Qscore etc.
            num_reads += 1
            read_sequence = str.upper(str(record.seq))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY_REV)
            if key_index >= 0: # if key found
                start_index = key_index + KEY_START
                guide = read_sequence[(start_index - 20):(start_index)]
                if guide in dict_perfects:
                    dict_perfects[guide] += 1
                    num_perfect_matches += 1
                else:
                    num_np_matches += 1
                    list_np.append(guide)
            else:
                num_nokey += 1
    else:
        raise Exception('ERROR! Specified direction is not valid')
    handle.close()

    # STEP 3: SORT DICTIONARIES AND GENERATE OUTPUT FILES
    # sort perf matches (A-Z) with guides,counts as k,v and output to csv
    df_perfects = pd.DataFrame(data=dict_perfects.items(), columns=['sgRNA_seq', 'reads'])
    df_perfects.sort_values(by='sgRNA_seq', inplace=True)
    df_perfects.to_csv(out_counts, index=False, header=False)
    # now sort non-perfect matches by frequency and output to csv
    dict_np = Counter(list_np) # use Counter to tally up np matches
    df_npmatches = pd.DataFrame(data=dict_np.items(), columns=['sgRNA_seq', 'reads'])
    df_npmatches.sort_values(by='reads', ascending=False, inplace=True)
    df_npmatches.to_csv(out_np, index=False)

    # STEP 4: CALCULATE STATS AND GENERATE STAT OUTPUT FILE
    # percentage of guides that matched perfectly
    pct_perfmatch = round(num_perfect_matches/float(num_perfect_matches + num_np_matches) * 100, 1)
    # percentage of undetected guides (no read counts)
    guides_with_reads = np.count_nonzero(list(dict_perfects.values()))
    guides_no_reads = len(dict_perfects) - guides_with_reads
    pct_no_reads = round(guides_no_reads/float(len(dict_perfects.values())) * 100, 1)
    # skew ratio of top 10% to bottom 10% of guide counts
    top_10 = np.percentile(list(dict_perfects.values()), 90)
    bottom_10 = np.percentile(list(dict_perfects.values()), 10)
    if top_10 != 0 and bottom_10 != 0:
        skew_ratio = top_10/bottom_10
    else:
        skew_ratio = 'Not enough perfect matches to determine skew ratio'
    # calculate the read coverage (reads processed / sgRNAs in library)
    num_guides = df_ref['sgRNA_seq'].shape[0]
    coverage = round(num_reads / num_guides, 1)
    # calculate the number of unmapped reads (num_nokey / total_reads)
    pct_unmapped = round((num_nokey / num_reads) * 100, 2)

    # write analysis statistics to statfile
    with open(out_stats, 'w') as statfile:
        statfile.write('Number of reads processed: ' + str(num_reads) + '\n')
        statfile.write('Number of reads where key was not found: ' + str(num_nokey) + '\n')
        statfile.write('Number of perfect guide matches: ' + str(num_perfect_matches) + '\n')
        statfile.write('Number of nonperfect guide matches: ' + str(num_np_matches) + '\n')
        statfile.write('Number of undetected guides: ' + str(guides_no_reads) + '\n')
        statfile.write('Percentage of unmapped reads (key not found): ' + str(pct_unmapped) + '\n')
        statfile.write('Percentage of guides that matched perfectly: ' + str(pct_perfmatch) + '\n')
        statfile.write('Percentage of undetected guides: ' + str(pct_no_reads) + '\n')
        statfile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio) + '\n')
        statfile.write('Read coverage: ' + str(coverage))
        statfile.close()

    print(str(in_fastq) + ' processed')
    return

#%% merge_and_norm() - merge sample read counts and log2/t0 normalize

def merge_and_norm(dict_counts, in_ref, t0='t0', dir_counts='', save='all',
                   out_folder='', out_reads='agg_reads.csv', out_log2='agg_log2.csv',
                   out_t0='agg_t0_reps.csv', return_df=None):
    """
    Aggregate raw counts and perform log2-transform and t0 normalization.

    For a given set of samples and their raw read count files from count_reads,
    aggregate them into a single dataframe, normalize counts to reads per million,
    perform log2 transform, and then normalize to t0. The aggregated raw reads,
    log2 transformed values, and t0 normalized values can be saved to csv.

    Parameters
    ----------
    dict_counts : dict in format {'sample name': 'file name'}
        Dictionary to map sample names (key) to read count file names (value),
        as key,value (e.g. {'KN-0': 'KN-0_counts.csv'}). Must include the t0
        sample in the dict.
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    t0 : str, default 't0'
        Name of the t0 sample in dict_counts. If you have multiple t0 samples
        (e.g. paired t0s for specific samples), then you will need to run this
        function separately for each set of samples with their appropriate t0.
    dir_counts : str, default ''
        Name of the subfolder to find the read count csv files. The default is
        the current working directory.
    save : {'all', None, ['reads', 'log2', 't0']}, default 'all'
        Choose files for export to csv. The default is 'all', which is the
        aggregated read counts ('reads'), log2 normalized values ('log2'), and
        t0 normalized values ('t0'). You may also enter any combination of
        'reads', 'log2', 't0') as a list of strings to choose which ones to
        save ('all' is equivalent to a list of all three). None will not export
        any files to csv.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_reads : str, default 'agg_reads.csv'
        Name of the aggregated raw reads csv output file.
    out_log2 : str, default 'agg_log2.csv'
        Name of the aggregated log2 normalized values csv output file.
    out_t0 : str, default 'agg_t0_reps.csv'
        Name of the aggregated t0 normalized values csv output file.
    return_df : {None, 'reads', 'log2', 't0'}, default None
        Whether to return a dataframe at function end. The default is None,
        which returns nothing. However, you can return the reads, log2 norm,
        or t0 norm values by calling 'reads', 'log2', or't0', respectively.
    """

    # import reference file, define variables, check for requirements
    path = Path.cwd()
    inpath = path / dir_counts
    df_ref = pd.read_csv(in_ref)
    if 'sgRNA_seq' not in df_ref.columns.tolist():
        raise Exception('in_ref is missing column: sgRNA_seq')
    if t0 not in dict_counts.keys():
        raise Exception ('dict_counts is missing the t0 sample')
    # rearrange the dict samples to place t0 first (for order purposes later)
    list_samples = [t0] + [samp for samp in dict_counts.keys() if samp != t0]
    # aggregate read counts from all samples into df_rawreads
    # also perform log2 norm (brian/broad method; log2(rpm + 1 / total reads))
    df_reads, df_log2, df_t0 = df_ref.copy(), df_ref.copy(), df_ref.copy()
    for sample in list_samples:
        df_temp = pd.read_csv(inpath / dict_counts[sample], names=['sgRNA_seq', sample])
        # aggregating raw reads
        df_reads = pd.merge(df_reads, df_temp, on='sgRNA_seq')
        # log2 normalization
        total_reads = df_reads[sample].sum()
        df_log2[sample] = df_reads[sample].apply(lambda x: np.log2((x * 1000000 / total_reads) + 1))
        # t0 normalization
        df_t0[sample] = df_log2[sample].sub(df_log2[t0])
    # drop the t0 column since it will be 0
    df_t0.drop(columns=t0, inplace=True)

    # export files and return dataframes if necessary
    outpath = path / out_folder
    Path.mkdir(outpath, exist_ok=True)
    # dictionary to map kws to dfs and output file names
    dict_df = {'reads': (df_reads, out_reads), 'log2': (df_log2, out_log2), 't0': (df_t0, out_t0)}
    # determine which files to export
    if save == 'all':
        save = ['reads','log2','t0']
    if isinstance(save, list):
        for key in save:
            dict_df[key][0].to_csv(outpath / dict_df[key][1], index=False)
    elif save is None:
        pass
    else:
        warnings.warn('Invalid value for save. No files exported')
    # determine df to return
    print('Merge and normalize completed')
    if return_df in dict_df.keys():
        return dict_df[return_df][0]
    elif return_df is None:
        return
    else:
        print('Invalid value for return_df. No dataframe returned')
        return

#%% average_reps() - average replicates by condition

def average_reps(dict_conds, in_lfc, in_ref, save=True, out_folder='',
                 out_conds='agg_t0_conds.csv', return_df=False):
    """
    For a set of conditions, average the replicates and export it to csv.

    Averages the replicates for each condition (e.g. treatment, control) and
    exports the csv file with averaged replicate values for each condition.
    Note that this function is generally used for averaging t0 normalized
    replicate values, but can be used on any set of replicate values
    (e.g. log2 fc values or even raw reads).

    Parameters
    ----------
    dict_conds : dict in format {'replicate': 'condition'}
        Dictionary to map replicates (key) to conditions (value). Replicates
        must match the column headers in the in_lfc file (i.e. sample names),
        otherwise they will be discarded before performing the averaging.
        Example: {'KN-1': 'unsorted', 'KN-2': 'unsorted', 'KN-3': 'sorted'}
    in_lfc : str or path
        String or path to the individual replicate values csv file. The column
        headers for the replicate values must match the keys in dict_conds
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    save : bool, default True
        Whether to save the averaged replicate values as a csv file.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_conds : str, default 'agg_t0_conds.csv'
        Name of the averaged replicate values csv output file.
    return_df : bool, default False
        Whether to return the averaged reps dataframe The default is False.
    """

    # import files, define variables, check for requirements, etc.
    path = Path.cwd()
    df_lfc = pd.read_csv(in_lfc)
    # make df to map replicates to condition
    df_map = pd.DataFrame(data=dict_conds.items(), columns=['rep','condition'])
    # check to make sure replicates are in the input file
    list_reps = df_map['rep'].tolist()
    if not all(rep in list_reps for rep in df_lfc.columns.tolist()):
        list_miss = [rep for rep in list_reps if rep not in df_lfc.columns.tolist()]
        # remove missing replicates from df_map and raise Warning
        df_map = df_map.loc[~df_map['rep'].isin(list_miss)]
        warnings.warn('in_lfc is missing replicates (removed from analysis): ' + str(list_miss))

    # generate df to hold the averaged replicates per condition
    df_conds = pd.read_csv(in_ref)
    for cond in df_map['condition'].unique().tolist():
        # for each unique condition, find the reps
        reps = df_map.loc[df_map['condition'] == cond]['rep'].tolist()
        # skip averaging for single replicates (otherwise breaks script)
        if len(reps) > 1:
            df_conds[cond] = df_lfc[reps].mean(axis=1)
        elif len(reps) == 1:
            df_conds[cond] = df_lfc[reps]
        else:
            raise Exception('Error! Replicate number not valid')

    # export files and return dataframes if necessary
    if save:
        outpath = path / out_folder
        Path.mkdir(outpath, exist_ok=True)
        df_conds.to_csv(outpath / out_conds, index=False)
    print('Average reps completed')
    if return_df:
        return df_conds
    else:
        return

#%% compare_conds() - compare treatment vs. control to calculate enrichments

def compare_conds(list_comparisons, in_lfc, in_ref, save=True, out_folder='',
                  out_comps='comparisons.csv', return_df=False):
    """
    Perform pairwise comparisons given a list and export the output to a csv.

    Given a list of comparisons (e.g. treatment vs. control), perform pairwise
    comparisons, generate a dataframe, and export to csv. The list of comparisons
    must be in the format (comparison name, condition 1, condition 2).
    The comparison is performed as (condition 1 - condition 2). Note that this
    can be applied to any format of values, not just averaged condition reps.

    Parameters
    ----------
    list_comparisons : list of tuples in format (name, sample 1, sample 2)
        A list of tuples denoting the comparisons to make, with the comparison
        being sample 1 - sample 2 (e.g. treatment - control). The output column
        headers will be labeled by the comparison name in the tuple.
    in_lfc : str or path
        String or path to the csv file containing the values for comparison.
        The column headers must match the sample names in list_comparisons
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    save : bool, default True
        Whether to save the comparisons dataframe as a csv file.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_comps : str, default 'comparisons.csv'
        Name of the comparisons csv output file.
    return_df : bool, default False
        Whether to return the comparisons dataframe. The default is False.
    """

    # import files, define variables, check for requirements
    path = Path.cwd()
    df_lfc = pd.read_csv(in_lfc)
    df_comps = pd.read_csv(in_ref)
    # perform treatment vs. control comparison
    for name, treatment, control in list_comparisons:
        df_comps[name] = df_lfc[treatment].sub(df_lfc[control])
    # export files and return dataframes if necessary
    if save:
        outpath = path / out_folder
        Path.mkdir(outpath, exist_ok=True)
        df_comps.to_csv(outpath / out_comps, index=False)
    print('Compare conditions completed')
    if return_df:
        return df_comps
    else:
        return

#%% batch_count() - batch process FASTQ files with input csv file

def batch_count(in_batch, in_ref, dir_fastq='', dir_counts='', dir_np='',
                dir_stats='', **kwargs):
    """
    Perform count_reads() on a batch scale given a sample sheet.

    Batch process FASTQ files with a csv sample sheet (in_batch) and a reference
    file (in_ref). in_batch must include headers for sample ids ('sample_id'),
    FASTQ file names ('fastq_file'), and treatment conditions ('condition').

    Parameters
    ----------
    in_batch : str or path
        String or path to the batch sample sheet csv file. Must have column
        headers for sample ids ('sample_id'), FASTQ file names ('fastq_file'),
        and treatment conditions (e.g. drug, DMSO; 'condition').
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    dir_fastq : str, default ''
        The subfolder containing the FASTQ files. The default is the current
        working directory.
    dir_counts : str, default ''
        The subfolder to export the sgRNA read count csv files. The default is
        the current working directory.
    dir_np : str, default ''
        The subfolder to export the non-perfect match csv files. The default is
        the current working directory.
    dir_stats : str, default ''
        The subfolder to export the read count statistics files. The default is
        the current working directory.
    **kwargs : key, value mappings in format x=y
        All other keyword arguments are passed to count_reads(). See the
        count_reads() documentation for more information. Other kwargs include:
        KEY_INTERVAL, DIR, KEY, KEY_REV.
    """

    batch_st = time.perf_counter()
    # define all the directory paths
    path = Path.cwd()
    list_dirs = [path / subdir for subdir in [dir_fastq, dir_counts, dir_np, dir_stats]]
    for subdir in list_dirs:
        Path.mkdir(subdir, exist_ok=True)

    # import batch csv and process samples with count_reads()
    df_batch = pd.read_csv(in_batch)
    list_reqcols = ['sample_id', 'fastq_file', 'condition']
    list_batchcols = df_batch.columns.tolist()
    if not all(col in list_batchcols for col in list_reqcols):
        list_miss = [col for col in list_reqcols if col not in list_batchcols]
        raise Exception('Error! in_batch is missing column(s): ' + str(list_miss))

    # perform batch processing
    for row in df_batch.itertuples():
        t_start = time.perf_counter()
        fastq = list_dirs[0] / row.fastq_file
        counts = list_dirs[1] / (row.sample_id + '_counts.csv')
        np = list_dirs[2] / (row.sample_id + '_npcounts.csv')
        stats = list_dirs[3] / (row.sample_id + '_stats.txt')
        count_reads(in_fastq=fastq, in_ref=in_ref, out_counts=counts,
                    out_np=np, out_stats=stats, **kwargs)
        t_end = time.perf_counter()
        print(row.sample_id + ' processed in %.2f sec' % (t_end - t_start))

    batch_end = time.perf_counter()
    print('Batch count completed in %.2f min' % ((batch_end - batch_st) / 60))
    return

#%% batch_process() - batch process count_reads output (merging, log2/t0 norm)

def batch_process(in_batch, in_ref, merge_stats=True, dir_counts='', dir_stats='',
                  in_counts=None, in_stats=None, save='all', out_folder='',
                  out_prefix='', return_df=None):
    """
    Batch merging and pre-processing (log2/t0) of samples given a sample sheet.

    Batch processing of read counts from count_reads using a csv sample sheet
    (in_batch) and a reference file (in_ref). Aggregates raw reads, performs
    log2 and t0 normalization, averages reps per condition for t0 normalized
    values, and exports to csv. Also merges and exports the stat files as csv.

    Parameters
    ----------
    in_batch : str or path
        String or path to the batch sample sheet csv file. Must have column
        headers for sample ids ('sample_id'), FASTQ file names ('fastq_file'),
        and treatment conditions (e.g. drug, DMSO; 'condition').
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    merge_stats : bool, default True
        Whether to merge the read counts statistics files.
    dir_counts : str, default ''
        The subfolder containing the read counts csv files. The default is
        the current working directory.
    dir_stats : str, default ''
        The subfolder containing the read counts stats files. The default is
        the current working directory.
    in_counts : list of tup in format ('sample_id', 'file name'), default None
        List of tuples of (sample_id, file name) for the samples
        in in_batch. The default is None, which assumes the default naming
        scheme from batch_count ('sample_id' + '_counts.csv' = 'KN-0_counts.csv').
        If your read counts files do not follow the default naming scheme,
        then use in_counts to map sample_ids to your read count file names.
    in_stats : list of tup in format ('sample_id', 'file name'), default None
        List of tuples of (sample_id, file name) for the samples
        in in_batch. The default is None, which assumes the default naming
        scheme from batch_count ('sample_id' + '_stats.txt' = 'KN-0_stats.txt').
        If your stats files do not follow the default naming scheme, then use
        in_stats to map sample_ids to your read count stats file names.
    save : {'all', None, ['reads', 'log2', 't0', 'conds', 'stats']}, default 'all'
        Choose files for export to csv. The default is 'all', which is the
        aggregated read counts ('reads'), log2 normalized values ('log2'), and
        t0 normalized values ('t0'). You may also enter any combination of
        'reads', 'log2', 't0') as a list of strings to choose which ones to
        save ('all' is equivalent to a list of all three). None will not export
        any files to csv.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_prefix : str, default ''
        File prefix for the output csv files. The prefix should contain an
        underscore.
    return_df : {None, 'reads', 'log2', 't0', 'conds', 'stats'}, default None
        Whether to return a dataframe at function end. The default is None,
        which returns nothing. However, you can return the reads, log2 norm,
        t0 norm, averaged reps by condition, or stats dataframes by calling
        'reads', 'log2', 't0', 'conds', or 'stats', respectively.
    """

    # import ref files and define variables/paths
    path = Path.cwd()
    df_ref = pd.read_csv(in_ref)
    if 'sgRNA_seq' not in df_ref.columns.tolist():
        raise Exception('in_ref is missing column: sgRNA_seq')
    df_batch = pd.read_csv(in_batch)
    list_reqcols = ['sample_id', 'fastq_file', 'condition']
    list_batchcols = df_batch.columns.tolist()
    if not all(col in list_batchcols for col in list_reqcols):
        list_miss = [col for col in list_reqcols if col not in list_batchcols]
        raise Exception('Error! in_batch is missing column(s): ' + str(list_miss))
    if 't0' not in df_batch['condition'].tolist():
        raise Exception('t0 condition not found in the in_batch file')
    # defaults to cwd if subdir == ''
    counts_path = path / dir_counts
    stats_path = path / dir_stats
    if in_counts is None:
        df_batch['counts_files'] = df_batch['sample_id'] + '_counts.csv'
    else:
        df_temp = pd.DataFrame(in_counts, columns=['sample_id', 'counts_files'])
        df_batch = df_batch.merge(df_temp, on='sample_id', how='left')
    if in_stats is None:
        df_batch['stats_files'] = df_batch['sample_id'] + '_stats.txt'
    else:
        df_temp = pd.DataFrame(in_stats, columns=['sample_id', 'stats_files'])
        df_batch = df_batch.merge(df_temp, on='sample_id', how='left')

    # import csv files and generate dfs for raw reads and log2 norm
    df_reads, df_log2 = df_ref.copy(), df_ref.copy()
    for row in df_batch.itertuples():
        file = counts_path / row.counts_files
        df_temp = pd.read_csv(file, names=['sgRNA_seq', row.sample_id])
        # merge on sgRNA_seq to aggregate columns
        df_reads = pd.merge(df_reads, df_temp, on='sgRNA_seq')
        # perform log2 normalization (brian/broad method)
        total_reads = df_reads[row.sample_id].sum()
        df_log2[row.sample_id] = df_reads[row.sample_id].apply(lambda x: np.log2((x * 1000000 / total_reads) + 1))

    # perform t0 normalization
    df_t0 = df_ref.copy()
    t0 = df_batch.loc[df_batch['condition'] == 't0']['sample_id']
    if t0.shape[0] != 1:
        raise Exception('Only a single t0 sample is allowed')
    t0 = t0[0]
    for row in df_batch.itertuples():
        df_t0[row.sample_id] = df_log2[row.sample_id].sub(df_log2[t0])
    df_t0.drop(columns=t0, inplace=True) # drop the t0 col

    # average replicates by condition
    list_conds = df_batch['condition'].unique().tolist()
    list_conds.remove('t0')
    df_conds = df_ref.copy()
    for cond in list_conds:
        reps = df_batch.loc[df_batch['condition'] == cond]['sample_id'].tolist()
        if len(reps) > 1:
            df_conds[cond] = df_t0[reps].mean(axis=1)
        elif len(reps) == 1:
            df_conds[cond] = df_t0[reps]
        else:
            raise Exception('Error! Invalid number of replicates')

    # merge statistics files
    if merge_stats:
        df_stats = pd.DataFrame(columns=['parameters'])
        for row in df_batch.itertuples():
            file = stats_path / row.stats_files
            df_temp = pd.read_csv(file, sep=': ', engine='python', names=['parameters', row.sample_id])
            df_stats = pd.merge(df_stats, df_temp, on='parameters', how='outer')

    # export files and return dataframes if necessary
    outpath = path / out_folder
    Path.mkdir(outpath, exist_ok=True)
    # dictionary to map kws to dfs and output file names
    dict_df = {'reads': (df_reads, out_prefix + 'reads.csv'),
               'log2': (df_log2, out_prefix + 'log2.csv'),
               't0': (df_t0, out_prefix + 't0_reps.csv'),
               'conds': (df_conds, out_prefix + 't0_conds.csv')}
    if merge_stats:
        dict_df.update({'stats': (df_stats, out_prefix + 'stats.csv')})
    # determine which files to export
    if save == 'all':
        save = ['reads','log2','t0', 'conds', 'stats']
    if isinstance(save, list):
        for key in save:
            dict_df[key][0].to_csv(outpath / dict_df[key][1], index=False)
    elif save is None:
        pass
    else:
        warnings.warn('Invalid value for save. No files exported')
    # determine df to return
    print('Batch processing completed')
    if return_df in dict_df.keys():
        return dict_df[return_df][0]
    elif return_df is None:
        return
    else:
        print('Invalid value for return_df. No dataframe returned')
        return

#%% merge_stats() - merge sample read count statistics files

def merge_stats(dict_stats, in_ref, dir_stats='', save=True, out_folder='',
                out_stats='agg_stats.csv', return_df=False):
    """
    Aggregates the stats files from count_reads() and outputs as a csv file.

    Parameters
    ----------
    dict_stats : dict in format {'sample name': 'file name'}
        Dictionary to map sample names (key) to stats files (value).
        For example, {'KN-0': 'KN-0_stats.txt'}
    dir_stats : str, default ''
        Name of the subfolder to find the stats text files. The default is the
        current working directory.
    save : bool, default True
        Whether to save the aggregated read count statistics as a csv file.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_stats : str, default 'agg_stats.csv'
        Name of the aggregated read count statistics csv output file.
    return_df : bool, default False
        Whether to return the aggregated stats dataframe. The default is False.
    """

    # import reference file, define variables, check for requirements
    path = Path.cwd()
    inpath = path / dir_stats
    # generate df for merging files
    df_stats = pd.DataFrame(columns=['parameters'])
    for sample, file in dict_stats.items():
        df_temp = pd.read_csv(inpath / file, sep=': ', engine='python', names=['parameters', sample])
        df_stats = pd.merge(df_stats, df_temp, on='parameters', how='outer')
    # transpose df_stats with samples as rows, params as columns
    df_stats = df_stats.transpose(copy=True)
    # export files and return dataframes if necessary
    if save:
        outpath = path / out_folder
        Path.mkdir(outpath, exist_ok=True)
        df_stats.to_csv(outpath / out_stats, index=True)
    print('Merge stats completed')
    if return_df:
        return df_stats
    else:
        return

#%% qc_stats() - rough QC of read count stats to see if they pass metrics

def qc_stats(dict_stats, dir_stats='', min_depth=100, target_depth=500,
             max_pct_unmapped=10, min_pct_perfect=70, max_pct_undetected=0.5,
             max_skew=10, out_folder='', save_report=True, save_df=False,
             out_report='QC_stats_report.txt', out_df='QC_stats.csv', return_df=False):
    """
    Check if the read count stats of a set of samples passes rough QC metrics

    Performs a rough QC check on the count_reads() output statistics files to 
    determine which samples do not pass the rough metrics for "good" quality.
    The metrics outlined in the Joung et al. Nature Protocols 2016 paper for a
    good library are: sequencing depth of 1000x, >70% perfect sgRNA matches,
    <0.5% undetected sgRNAs, and a skew ratio <10. On top of this, I've added
    a metric for max % of unmapped reads (key not found; default 10%). The
    default rough QC metrics generally follow the Joung guidelines for min % of
    perfect matches (default: 70%), max % of undetected sgRNAs (default: 0.5%),
    and max skew ratio (default: 10), but I've lowered the default target read
    depth to 500x and added a minimum depth of 100x (b/w 100x-500x is "OK").

    Nota bene: samples 'failing' these metrics are not necessarily BAD samples,
    and in fact, samples that fail these metrics generally suggest that a screen
    is working (e.g. late screen timepoints should have highly asymmetric sgRNA
    distributions (large skew), sgRNAs that drop out (large # of undetected), etc.)
    The QC stats report gives a rough overview of NGS and sample quality.

    Parameters
    ----------
    dict_stats : dict in format {'sample name': 'file name'}
        Dict to map sample names to stats files, with sample name, file name
        as key,value (e.g. {'KN-0': 'KN-0_stats.txt'}).
    dir_stats : str, default ''
        Name of the subfolder to find the stats text files. The default is the
        current working directory.
    min_depth : int, default 100
        The minimum read depth that a sample should exceed. This is calculated
        as the # of total reads / # of sgRNAs in the library. The default is
        100x, but can be as low as 50x.
    target_depth : int, default 500
        The target read depth that a sample should meet. This is calculated
        as the # of total reads / # of sgRNAs in the library. The default is
        500x, but can be lower. Samples that fall within the min_depth and the
        target_depth are considered as "OK" and having moderate read coverage.
    max_pct_unmapped : int or float, default 10
        The max percentage of total reads that are unmapped (no key found).
        This is calculated as # of nokey reads / # of total reads. The default
        is 10% (as my samples generally range between 2-6%). Samples exceeding
        this threshold likely have an issue with the count_reads() parameters
        or have poor NGS read quality.
    min_pct_perfect : int or float, default 70
        The minimum percentage of mapped reads that match an sgRNA in the in_ref
        file. Calculated as the # of perf matches / # reads with key. The
        default is 70%.
    max_pct_undetected : int or float, default 0.5
        The max percentage of sgRNAs in the library with no reads. Calculated
        as the # of sgRNAs with zero reads / # of total sgRNAs in library. The
        default is 0.5%, but can be adjusted such as when the library is small
        and a single missing sgRNA may exceed the 0.5% threshold.
    max_skew : int, 10
        The max skew ratio allowed. The skew ratio (not true skew) is defined
        as the # reads in the 90th percentile / # reads in the 10th percentile.
        The default is 10. Note that uncalculated skew ratios ("Not enough
        perfect matches...") always fail the max_skew test.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    save_report : bool, default True
        Whether to save the QC stats report. The default is True.
    save_df : bool, default False
        Whether to save the QC stats dataframe as a csv. The default is False.
    out_report : str, default 'QC_stats_report.txt'
        The name of the output QC stats report .txt file.
    out_df : str, default 'QC_stats.csv'
        The name of the output QC stats dataframe .csv file
    return_df : bool, default False
        Whether to return the final QC dataframe. The default is False.
    """

    # import reference file, define variables, check for requirements
    path = Path.cwd()
    inpath = path / dir_stats
    list_params = [('Number of reads processed', 'total_reads'),
                   ('Number of reads where key was not found', 'n_unmapped'),
                   ('Number of perfect guide matches', 'n_perfect'),
                   ('Number of nonperfect guide matches', 'n_nonperfect'),
                   ('Number of undetected guides', 'n_undetected'),
                   ('Percentage of unmapped reads (key not found)', 'pct_unmapped'),
                   ('Percentage of guides that matched perfectly', 'pct_perfect'),
                   ('Percentage of undetected guides', 'pct_undetected'),
                   ('Skew ratio of top 10% to bottom 10%', 'skew_ratio'),
                   ('Read coverage', 'coverage')]
    # iterate through stats files and aggregate into df_stats
    df_stats = pd.DataFrame(data=list_params, columns=['parameters','abbr'])
    for sample, file in dict_stats.items():
        df_temp = pd.read_csv(inpath / file, sep=': ', engine='python', names=['parameters', sample])
        df_stats = pd.merge(df_stats, df_temp, on='parameters', how='outer')
    df_stats.set_index(keys=df_stats['abbr'], inplace=True)
    df_stats.drop(columns=['parameters', 'abbr'], inplace=True)
    # convert dataframe into numeric dtypes for processing
    df_stats = df_stats.transpose(copy=True)
    list_dtypes = [int,int,int,int,int,float,float,float,float,float]
    dict_dtypes = {col:dtype for col,dtype in zip(df_stats.columns.tolist(), list_dtypes)}
    df_stats = df_stats.astype(dtype=dict_dtypes, errors='ignore')
    df_stats['skew_ratio'] = df_stats['skew_ratio'].apply(lambda x: np.where(isinstance(x, str), 0, x))
    df_stats['skew_ratio'] = df_stats['skew_ratio'].astype(float)
    # generate df to check for QC metrics
    df_qc = pd.DataFrame(index=df_stats.index)
    # read coverage - pass if >= ideal_depth, OK if b/w min and ideal, fail if < min_depth
    df_qc['coverage_qc'] = np.where(df_stats['coverage'] >= target_depth, 'pass',
                                    (np.where(df_stats['coverage'] < min_depth, 'fail', 'ok')))
    # unmapped (key not found) reads - pass if < max_pct_unmapped (default: 10%)
    df_qc['unmapped_qc'] = np.where(df_stats['pct_unmapped'] < max_pct_unmapped, 'pass', 'fail')
    # perfect matches - pass if > min_pct_perfect (default: 70%)
    df_qc['perfects_qc'] = np.where(df_stats['pct_perfect'] > min_pct_perfect, 'pass', 'fail')
    # undetected sgRNAs - pass if < pct_undetected (default: 0.5%)
    df_qc['undetected_qc'] = np.where(df_stats['pct_undetected'] < max_pct_undetected, 'pass', 'fail')
    # skew ratio - pass if < max_skew (default: 10)
    # 0 is stand-in for NaN (not enough matches), which is fail
    df_qc['skew_qc'] = np.where(df_stats['skew_ratio'] > max_skew, 'fail',
                                (np.where(df_stats['skew_ratio'] != 0, 'pass', 'fail')))

    # generate fastqc report and export files
    outpath = path / out_folder
    Path.mkdir(outpath, exist_ok=True)
    if save_report:
        coverage = '(' + str(min_depth) + 'x-' + str(target_depth) + 'x)'
        list_modcov = df_qc.loc[df_qc['coverage_qc'] == 'ok'].index.tolist()
        list_failcov = df_qc.loc[df_qc['coverage_qc'] == 'fail'].index.tolist()
        list_unmapped = df_qc.loc[df_qc['unmapped_qc'] == 'fail'].index.tolist()
        list_perfects = df_qc.loc[df_qc['perfects_qc'] == 'fail'].index.tolist()
        list_undetected = df_qc.loc[df_qc['undetected_qc'] == 'fail'].index.tolist()
        list_skew = df_qc.loc[df_qc['skew_qc'] == 'fail'].index.tolist()
        with open(outpath / out_report, 'w') as outfile:
            outfile.write('Samples with moderate read coverage ' + coverage + ': ' + str(list_modcov) + '\n')
            outfile.write('Samples with <' + str(min_depth) + 'x read coverage: ' + str(list_failcov) + '\n')
            outfile.write('Samples with >' + str(max_pct_unmapped) + '% unmapped reads: ' + str(list_unmapped) + '\n')
            outfile.write('Samples with <' + str(min_pct_perfect) + '% perfect matches: ' + str(list_perfects) + '\n')
            outfile.write('Samples with >' + str(max_pct_undetected) + '% undetected sgRNAs: ' + str(list_undetected) + '\n')
            outfile.write('Samples with uncalculated or >' + str(max_skew) + ' skew ratio: ' + str(list_skew))
            outfile.close()
    if save_df:
        df_qc.to_csv(outpath / out_df)
    print('QC stats completed')
    if return_df:
        return df_qc
    else:
        return

#%% qc_plot_dist() - plot sample read count distribution as hist and Lorenz curve

def qc_plot_dist(sample, in_count=None, scale='talk', title=None, save_as='pdf',
                 out_folder='', out_fig='dist_plot', return_fig=True):
    """
    Analyze a sample's sgRNA distribution and plot a histogram and Lorenz curve.

    QC function to analyze the sgRNA distribution of a sample (generally the
    sgRNA library) and visualizes it as a histogram of the read counts as well
    as a Lorenz curve (P-P plot of cdf(reads) vs cdf(sgRNA, abundance-ranked)).
    Takes the count_reads perfect matches csv file ('counts.csv') as input. The
    Lorenz curve generally gives an idea of how even the distribution is based
    on how close it matches the ideal (straight line). Gini coefficients can be
    calculated from the AUC values given in the graph, but usually unnecessary.

    Parameters
    ----------
    sample : str
        The name of the sample to be analyzed. This is used for labeling and to
        find the read counts csv file if no input is given (see in_count)
    in_count : str or path, default None
        String or path to the read counts csv file. If not specified, it defaults
        to sample + '_counts.csv' (e.g. 'KN-0_counts.csv') using the user's 
        sample name input and assumes the file is in the current working directory.
    scale : {'paper', 'notebook', 'talk', 'poster'}, default 'talk'
        The scaling factor to be used. This input is passed to sns.set_context().
        Scales the default mpl sizing by {0.8x, 1.0x, 1.3x, 1.6x} for {paper,
        notebook, talk, poster}, respectively. The default is 'talk' (1.3x).
    title : str, default None
        The title of the figure. The default is None, which will default to
        sample + ' sgRNA Distribution' (e.g. KN-0 sgRNA Distribution).
    save_as : {['pdf', 'png', etc.], None}, default 'pdf'
        String of the format to save the figure. Saving the figure in multiple
        formats is done by passing a list of strings of the desired formats
        (e.g. pass ['pdf','png'] to save a .pdf and a .png of the figure). The
        default is 'pdf', but any format supported by mpl.savefig() is allowed.
        If save_as is None, will not save any files.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_fig : str, default 'dist_plot'
        Name of the output figure file(s). NOTE: out_fig should not have a file
        extension as that is determined by the save_as input. For example, the
        default is 'dist_plot', which becomes 'dist_plot.pdf' if save_as='pdf'.
    return_fig : bool, default True
        Whether to return the generated figure. The default is True.
    """
    # if no input is specified for in_count, use default naming scheme
    if in_count is None:
        in_count = sample + '_counts.csv'
    # import read counts csv file
    df_dist = pd.read_csv(in_count, names=['sgRNA_seq', sample])
    # calculate cumulative frequency (Lorenz curve)
    df_dist = df_dist.sort_values(by=[sample], ascending=False).reset_index(drop=True)
    df_dist['read_fraction'] = df_dist[sample] / df_dist[sample].sum()
    df_dist['cum_sgRNAs'] = [i+1 for i in range(0, df_dist.shape[0])]
    df_dist['cum_reads'] = df_dist['read_fraction'].cumsum()
    # calculate the area under the curve (AUC) with np.trapz
    AUC = round(np.trapz(df_dist['cum_reads']) / df_dist.shape[0], 2)

    # set up style parameters for plotting in seaborn/mpl
    sns.set_context(scale)
    # plot histogram of sgRNA distribution
    fig, (ax1,ax2) = plt.subplots(ncols=2, figsize=(16,9))
    sns.distplot(df_dist[sample], ax=ax1, bins=50, kde=False, color='#c44e52', hist_kws={'alpha': 1})
    # plot the lorenz curve (cdf(reads) vs. cdf(sgRNAs) aka P-P plot)
    sns.lineplot(x='cum_sgRNAs', y='cum_reads', data=df_dist, ax=ax2, color='#c44e52',
                 label=sample + ' (AUC = ' + str(AUC) + ')')
    # plot the ideal curve (line of equality; y step is same for all x)
    ax2.plot([0, df_dist.shape[0]], [0,1], c='k', ls='--', label='Ideal (AUC = 0.50)')
    # aesthetics (histogram)
    ax1.set(xlim=(0, None), ylim=(0, None), xlabel='# of reads', ylabel='# of sgRNAs')
    # aesthetics (Lorenz curve)
    ax2.set(xlim=(0, df_dist.shape[0]), ylim=(0, 1), xlabel='sgRNAs ranked by abundance',
            ylabel='Cumulative fraction of reads')
    ax2.legend(loc='lower right')
    # aesthetics (global)
    for ax in ax1,ax2:
        for edge,spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_color('k')
    if title is None:
        title = sample + ' sgRNA Distribution'
    fig.suptitle(title)

    # exporting and returning plots
    outpath = Path.cwd() / out_folder
    Path.mkdir(outpath, exist_ok=True)
    if isinstance(save_as, str): # convert save_as to list
        save_as = [save_as]
    if isinstance(save_as, list):
        for ftype in save_as:
            out_plot = outpath / (out_fig + '.' + ftype)
            fig.savefig(out_plot, format=ftype, bbox_inches='tight', dpi=300)
    elif save_as is None:
        pass
    else:
        warnings.warn('Invalid value for save_as. Figure not exported.')
    print('sgRNA distribution plotting completed')
    if return_fig:
        return fig
    else:
        plt.close()
        return

#%% qc_pearson() - calculate Pearson correlations for samples and plot heatmap

def qc_pearson(list_samples, in_lfc, scale='talk', size=(12,9), title='Pairwise correlation heatmap',
               cmap='KCN_Reds', vlims=(-1,1), save_df=True, savefig_as='pdf',
               out_folder='', out_file='QC_pearson', return_obj='fig',
               plot_kws=dict(linewidths=0.5, linecolor='k', square=True)):
    """
    Calculate the Pearson correlation coefficient matrix and plot a heatmap.

    QC function to calculate the Pearson correlation coefficient matrix for a
    given list of samples (list_samples) and enrichment values (in_lfc) and
    generates a heatmap visualization of the pairwise correlation matrix. The
    list of samples must match the column headers in the input values file. This
    function is generally used to analyze the consistency between replicates
    (biological/technical) using the aggregated log2-transformed read counts csv
    output from merge_and_norm() or batch_process(), but any set of samples/lfc
    values can be used (for example, using averaged reps to compare conditions
    if you have multiple drug treatments).

    Parameters
    ----------
    list_samples : list of strings
        List of sample names to calculate the Pearson correlation coefficient
        matrix. Must match the column headers in the input file (in_lfc).
    in_lfc : str or path
        String or path to the csv file containing the enrichment values. The
        column headers must match the sample names in list_samples
    scale : {'paper', 'notebook', 'talk', 'poster'}, default 'talk'
        The scaling factor to be used. This is passed to sns.set_context().
        Scales the default mpl sizing by {0.8x, 1.0x, 1.3x, 1.6x} for {paper,
        notebook, talk, poster}, respectively. The default is 'talk' (1.3x).
    size : tuple of ints in format (width, height), default (12,9)
        The size of the figure in inches in the format (width, height). This is
        passed to plt.subplots. By default, the figure is made in 4:3 format,
        but it may be useful to use different dimensions (e.g. 16:9) for
        different applications. However, this may affect the sizing ratio of
        the objects being plotted.
        lowest color in the colormap.
    title : str, default 'Pairwise correlation heatmap'
        The title of the figure.
    cmap : str, default 'KCN_Reds'
        The color palette for the heatmap. This is passed to sns.heatmap(), which
        accepts matplotlib colormap names/objects or a custom palette in the mpl
        format. The default is 'KCN_Reds', which is a custom monocolor gradient
        starting from grey96 to red. The gradient is similar to the built-in
        matplotlib cmap 'Reds', except the base is grey rather than light red.
    vlims : tuple of floats in format (vmin, vmax), default (-1,1)
        The min and max values for the heatmap scale. The default is -1 to 1,
        which is the maximum range of the Pearson correlation coefficient.
        Setting a min/max smaller than the full range of values will cause all
        values exceeding that threshold to be colored as the min/max. For
        example, setting a min at 0 will cause all values <= 0 to map to the
    save_df : bool, default True
        Whether to save the Pearson correlation coefficient matrix.
    savefig_as : {['pdf', 'png', etc.], None}, default 'pdf'
        String of the format to save the figure. Saving the figure in multiple
        formats is done by passing a list of strings of the desired formats
        (e.g. pass ['pdf','png'] to save a .pdf and a .png of the figure). The
        default is 'pdf', but any format supported by mpl.savefig() is allowed.
        If save_figas is None, it will not save the figure.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_file : str, default 'QC_pearson'
        Name of the output file(s) for both the matrix (df) and the figures.
        NOTE: out_file should not have a file extension. For the matrix, the
        output will be csv, but the figures will have varying output types based
        on the format specified in savefig_as. For example, the default is
        'QC_pearson', which will return 'QC_pearson.csv' for the matrix and
        'QC_pearson.pdf' if savefig_as='pdf'.
    return_obj : {'fig', 'df', None}
        Whether to return an object generated by the function. The options are
        the Pearson matrix ('df'), the figure ('fig'), or nothing (None).
        The default is 'fig'.
    plot_kws : dict of {kwarg: value} mappings
        Keyword arguments passed to sns.heatmap() to adjust plot properties that
        are not covered in the preceding arguments. These kwargs must be
        interpretable by sns.heatmap(). By default, the kwargs passed in
        plot_kws are linewidths=0.5, linecolor='k', and square=True. Linewidths
        and linecolor set the width and color of the lines dividing the cells,
        respectively, and square=True sets the shape of the heatmap to make each
        cell a square. See the sns.heatmap() and the sns/mpl documentation for
        more information. NOTE: kwargs passed through plot_kws should not be
        previously defined. For example, do not pass vmin/vmax or cmap through
        plot_kws as these are set by vlims and cmap already.
    """

    # import lfc values from the input file (in_lfc)
    df_lfc = pd.read_csv(in_lfc)
    # generate df for pairwise Pearson correlation coefficient matrix
    df_pearson = pd.DataFrame(data=np.corrcoef(x=df_lfc[list_samples], rowvar=False),
                              columns=list_samples, index=list_samples)
    df_pearson = df_pearson.round(decimals=3)

    # set up style parameters for plotting in seaborn/mpl
    sns.set_context(scale)
    # if cmap == 'KCN_Reds', generate the mpl colormap object
    if cmap == 'KCN_Reds':
        # grey96 + the mpl 'Reds' colormap colors
        colors = ['#f5f5f5'] + list(sns.color_palette('Reds').as_hex())
        cmap = mpl.colors.LinearSegmentedColormap.from_list('KCN_Reds', colors=colors)
    # generate heatmap w/ sns.heatmap
    fig, ax = plt.subplots(figsize=size)
    sns.heatmap(data=df_pearson, ax=ax, cmap=cmap, vmin=vlims[0], vmax=vlims[1], **plot_kws)
    # aesthetics
    fig.suptitle(title)

    # exporting and returning files
    outpath = Path.cwd() / out_folder
    Path.mkdir(outpath, exist_ok=True)
    if save_df:
        df_pearson.to_csv(outpath / (out_file + '.csv'))
    if isinstance(savefig_as, str): # convert savefig_as to list
        savefig_as = [savefig_as]
    if isinstance(savefig_as, list):
        for ftype in savefig_as:
            out_plot = outpath / (out_file + '.' + ftype)
            fig.savefig(out_plot, format=ftype, bbox_inches='tight', dpi=300)
    elif savefig_as is None:
        pass
    else:
        warnings.warn('Invalid value for savefig_as. Figure not exported.')
    print('Pearson correlation calculation and plotting completed')
    dict_obj = {'fig': fig, 'df': df_pearson}
    if isinstance(return_obj, str):
        return dict_obj[return_obj]
    elif return_obj is None:
        plt.close()
        return
    else:
        print('Invalid value for return_obj. No object returned.')
        plt.close()
        return

#%% qc_pca() - perform PCA on samples and plot as scatter plot

def qc_pca(list_samples, in_lfc, dict_groups=None, group_label=None, scale='talk',
           colors='deep', size=(16,9), title='PCA scatter plot', save_as='pdf',
           out_folder='', out_fig='dist_plot', return_fig=True, **kwargs):
    """
    Perform PCA on a set of samples and plot it as a scatter plot.

    QC function to perform PCA (2 dimensions) on a given set of samples
    (e.g. replicates) and visualizes the resulting data as a scatter plot.
    Generally, this is used to visualize the similarity between samples and
    to see (or confirm) which samples are more closely correlated than others.
    For example, performing this as a QC check on a list of samples before
    averaging should show tight clustering between technical/biological
    replicates. Although this function is designed to check the correlations of
    technical replicates, it can be applied to any set of samples and enrichment
    values.

    Parameters
    ----------
    list_samples : list of strings
        List of sample names to calculate the perform PCA on. The sample names
        must match the column headers in the input file (in_lfc).
    in_lfc : str or path
        String or path to the csv file containing the enrichment values. The
        column headers must match the sample names in list_samples.
    dict_groups : dict of {'sample name': 'group'}, default None
        Dict to map samples (key) to 'groups' (value) to color the scatter plot
        points. A group is any classifier upon which samples can be grouped
        into, such as treatment, time, etc. For example, 3 technical replicates
        of +drug at t=21d could be grouped under 'drug_21d' as follows:
        {'rep1': 'drug_21d', 'rep2': 'drug_21d', 'rep3': 'drug_21d'}. The
        default is None, which colors each sample individually, but this is not
        recommended for large data sets.
    group_label : str, default None
        The name of the 'group' used to categorize samples in dict_groups. The
        group_label is also used to label the legend. The function will raise
        an Exception if a dict_groups input is provided, but no group_label is
        defined. The default is None, which sets the label to 'Sample.'
    scale : {'paper', 'notebook', 'talk', 'poster'}, default 'talk'
        The scaling factor to be used. This is passed to sns.set_context().
        Scales the default mpl sizing by {0.8x, 1.0x, 1.3x, 1.6x} for {paper,
        notebook, talk, poster}, respectively. The default is 'talk' (1.3x).
    colors : str of a palette name, None, list, or dict
        The color palette for the scatter plot. If a string is passed to colors,
        it must be interpretable by mpl/seaborn (see sns.color_palette()).
        The default is 'deep', which is the seaborn deep palette. If None is
        passed, then the function will use my custom mpl color palette.
    size : tuple of ints in format (width, height), default (16,9)
        The size of the figure in inches in the format (width, height). This is
        passed to plt.subplots. By default, the figure is made in 16:9 widescreen
        format, but it may be useful to use different dimensions (e.g. 4:3) for
        different applications. However, this may affect the sizing ratio of
        the objects being plotted.
    title : str, default 'PCA scatter plot'
        The title of the figure.
    save_as : {['pdf', 'png', etc.], None}, default 'pdf'
        String of the format to save the figure. Saving the figure in multiple
        formats is done by passing a list of strings of the desired formats
        (e.g. pass ['pdf','png'] to save a .pdf and a .png of the figure). The
        default is 'pdf', but any format supported by mpl.savefig() is allowed.
        If save_as is None, will not save any files.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_fig : str, default 'dist_plot'
        Name of the output figure file(s). NOTE: out_fig should not have a file
        extension as that is determined by the save_as input. For example, the
        default is 'dist_plot', which becomes 'dist_plot.pdf' if save_as='pdf'.
    return_fig : bool, default True
        Whether to return the generated figure. The default is True.
    **kwargs : key, value mappings in format x=y
        Other keyword arguments are passed to sns.scatterplot(). See the seaborn
        documentation for more information.
    """

    # import lfc values from the input file (in_lfc)
    df_lfc = pd.read_csv(in_lfc)
    # isolate values and transpose (samples as rows, sgRNAs as cols)
    df_lfc = df_lfc.loc[:, list_samples].transpose()
    # perform feature scaling (z-score normalization): ((x - mean) / std)
    df_scaled = StandardScaler().fit_transform(df_lfc)
    # perform PCA to reduce df_scaled to 2 dimensions and turn into df
    pca = PCA(n_components=2)
    df_pca = pca.fit_transform(df_scaled)
    # explained variance for PC1 and PC2
    var_PC1 = round(pca.explained_variance_ratio_[0] * 100, 2)
    var_PC2 = round(pca.explained_variance_ratio_[1] * 100, 2)
    # convert the np array to df with sample names as index
    df_pca = pd.DataFrame(data=df_pca, columns=['PC1','PC2'])
    df_pca['samples'] = list_samples

    # add col to color/label legend on dict_groups; if None, color by sample
    if group_label is None:
        if dict_groups is None:
            group_label = 'Sample'
        else: # if dict_groups supplied but no group_label
            raise Exception('A group_label is required if dict_groups is provided')
    if dict_groups is None:
        df_pca[group_label] = list_samples
    elif isinstance(dict_groups, dict):
        df_map = pd.DataFrame(data=dict_groups.items(), columns=['samples', group_label])
        df_pca = df_pca.merge(df_map, on='samples', how='left')
    else:
        raise Exception('Invalid value for dict_groups.')
    # if colors is None, use my custom mpl palette (default is seaborn deep)
    # it's better than the actual mpl default palette, trust me
    if colors is None:
        colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#999999',
                  '#fb9a99', '#a6cee3', '#b2df8a', '#cab2d6', '#fdbf6f', '#ffff33']
                  # red, blue, green, purple, orange, grey,
                  # lt red, lt blue, lt green, lt purple, lt orange, yellow
    # check if there are a sufficient number of unique colors, else raise Warning
    if isinstance(colors, str):
        num_colors = len(sns.color_palette(colors))
    else:
        num_colors = len(colors)
    if num_colors < len(df_pca[group_label].unique().tolist()):
        warnings.warn('There are more groups than colors; some groups will share colors.')

    # set up style parameters for plotting in seaborn/mpl
    sns.set_context(scale)
    # visualize PCA on scatterplot
    fig, ax = plt.subplots(figsize=size)
    sns.scatterplot(ax=ax, data=df_pca, x='PC1', y='PC2', hue=group_label,
                    palette=colors, **kwargs)
    # aesthetics
    ax.set(xlabel='PC1 (' + str(var_PC1) + '%)', ylabel='PC2 (' + str(var_PC2) + '%)')
    ax.grid(b=True, which='major', axis='both', lw=0.5, c='k', alpha=0.2)
    for edge,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    fig.suptitle(title)
    # by default mpl puts legend inside plot; set legend outside plot center-right
    # center left point of the legend box set to x=1.02, y=0.5 (plot is [0,1] for xy)
    ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))

    # exporting and returning files
    outpath = Path.cwd() / out_folder
    Path.mkdir(outpath, exist_ok=True)
    if isinstance(save_as, str): # convert save_as to list
        save_as = [save_as]
    if isinstance(save_as, list):
        for ftype in save_as:
            out_plot = outpath / (out_fig + '.' + ftype)
            fig.savefig(out_plot, format=ftype, bbox_inches='tight', dpi=300)
    elif save_as is None:
        pass
    else:
        warnings.warn('Invalid value for save_as. Figure not exported.')
    print('PCA plot completed')
    if return_fig:
        return fig
    else:
        plt.close()
        return