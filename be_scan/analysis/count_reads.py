"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Count the reads in a FASTQ file and assign them to a reference sgRNA set}
"""

from collections import Counter
import gzip
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

def count_reads(sample_sheet, in_ref, 
                file_dir='', 
                KEY_INTERVAL=(10,80), KEY='CGAAACACC', KEY_REV='GTTTTAGA', 
                dont_trim_G=False,
                out_dir='', out_file='counts_library.csv',
                save=True, return_df=True, save_files=True,
                ):
    """
    Given a set of sgRNA sequences and a FASTQ file, count the reads in the
    FASTQ, assign the reads to sgRNAs, and export the counts to a csv file `out_counts`. All
    sgRNA sequences not found in the reference file (non-perfect matches) are
    written to a separate csv file `out_nc`.

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'fastq_file', 'counts_file', 'noncounts_file', 'stats_file'
        a sheet with information on sequence id, 
        in_fastq (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_nc (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.

    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_REV : str, default 'GTTTTAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.

    file_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'counts_library.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe
    save_files : bool, default True
        Whether or not to save individual counts, noncounts, and stats files
    """
    sample_sheet = Path(sample_sheet)
    df = pd.read_csv(sample_sheet)
    samples = [list(a) for a in zip(df.fastq_file, df.counts_file, df.noncounts_file, df.stats_file)]

    # STEP 1A: OPEN INPUT FILES FOR PROCESSING, CHECK FOR REQUIRED FORMATTING
    # look for 'sgRNA_seq' column, raise Exception if missing
    in_ref = Path(in_ref)
    df_ref = pd.read_csv(in_ref, header=0) # explicit header = first row
    if 'sgRNA_seq' not in df_ref.columns.tolist():
        raise Exception('in_ref is missing column: sgRNA_seq')

    for fastq, counts, nc, stats in samples: 
        # fastq file of reads and paths to all output files, imported from sample_sheet
        in_fastq = Path(file_dir) / fastq
        out_counts, out_nc, out_stats = Path(out_dir) / counts, Path(out_dir) / nc, Path(out_dir) / stats        
        # try opening input FASTQ, raise Exception if not possible
        handle = gzip.open(in_fastq, 'rt') if str(in_fastq).endswith('.gz') else open(in_fastq, 'rt')

        # STEP 1B: SET UP VARIABLES FOR SCRIPT
        # make dictionary to hold sgRNA counts - sgRNA_seq, count as k,v
        dict_p = {sgRNA:0 for sgRNA in df_ref['sgRNA_seq']}
        list_np = [] # placeholder list for non-perfect matches
        # reads count of: total, perfect match, non perfect match, no key found, not 20bps
        num_reads, num_p_matches, num_np_matches, num_nokey, num_badlength = 0, 0, 0, 0, 0
        KEY_START, KEY_END = KEY_INTERVAL[0], KEY_INTERVAL[1] # set the key interval

        # STEP 2: PROCESS FASTQ FILE READS AND ADD COUNTS TO DICT
        while True: # contains the seq and Qscore etc.
            read = handle.readline()
            if not read: # end of file
                break
            elif read.startswith("@"): # if line is a read header
                read = handle.readline() # next line is the actual sequence
            else:
                continue
            num_reads += 1
            read_sequence = str.upper(str(read))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY)
            key_rev_index = key_region.find(KEY_REV)
            if key_index < 0 or key_rev_index <= key_index: # if keys not found
                num_nokey += 1
                continue
            start_index = key_index + KEY_START + len(KEY)
            end_index = key_rev_index + KEY_START
            guide = read_sequence[start_index:end_index]
            if not dont_trim_G:
                if guide.startswith('G') and len(guide) == 21:
                    guide = guide[1:]
            if len(guide) != 20:
                num_badlength += 1
                continue
            if guide in dict_p:
                dict_p[guide] += 1
                num_p_matches += 1
            else:
                num_np_matches += 1
                list_np.append(guide)
        handle.close()

        # STEP 3: SORT DICTIONARIES AND GENERATE OUTPUT FILES
        # sort perf matches (A-Z) with guides, counts as k,v and output to csv
        counts_name = counts.split("/")[-1]
        df_perfects = pd.DataFrame(data=dict_p.items(), columns=['sgRNA_seq', counts_name])
        if save_files:
            df_perfects.sort_values(by=counts_name, ascending=False, inplace=True)
            df_perfects.to_csv(out_counts, index=False)
        # add matching counts to dataframe
        df_ref = pd.merge(df_ref, df_perfects, on='sgRNA_seq', how='outer')
        df_ref[counts_name] = df_ref[counts_name].fillna(0)
        # now sort non-perfect matches by frequency and output to csv
        dict_np = Counter(list_np) # use Counter to tally up np matches
        nc_name = nc.split("/")[-1]
        df_ncmatches = pd.DataFrame(data=dict_np.items(), columns=['sgRNA_seq', nc_name])
        if save_files:
            df_ncmatches.sort_values(by=nc_name, ascending=False, inplace=True)
            df_ncmatches.to_csv(out_nc, index=False)
        # calculate the read coverage (reads processed / sgRNAs in library)
        num_guides = df_ref['sgRNA_seq'].shape[0]
    
        # STEP 4: CALCULATE STATS AND GENERATE STAT OUTPUT FILE
        # percentage of guides that matched perfectly
        pct_p_match = round(num_p_matches/float(num_p_matches + num_np_matches) * 100, 1)
        # percentage of undetected guides (no read counts)
        guides_with_reads = np.count_nonzero(dict_p.values())
        guides_no_reads = len(dict_p) - guides_with_reads
        pct_no_reads = round(guides_no_reads/float(len(dict_p.values())) * 100, 1)
        # skew ratio of top 10% to bottom 10% of guide counts
        top_10 = np.percentile(list(dict_p.values()), 90)
        bottom_10 = np.percentile(list(dict_p.values()), 10)
        if top_10 != 0 and bottom_10 != 0:
            skew_ratio = top_10/bottom_10
        else:
            skew_ratio = 'Not enough perfect matches to determine skew ratio'
        # calculate the read coverage (reads processed / sgRNAs in library)
        coverage = round(num_reads / num_guides, 1)
        # calculate the number of unmapped reads (num_nokey / total_reads)
        pct_unmapped = round((num_nokey / num_reads) * 100, 2)
        # write analysis statistics to statfile
        if save_files:
            with open(out_stats, 'w') as statfile:
                statfile.write('Number of reads processed: ' + str(num_reads) + '\n')
                statfile.write('Number of reads where key was not found: ' + str(num_nokey) + '\n')
                statfile.write('Number of reads where length was not 20bp: ' + str(num_badlength) + '\n')
                statfile.write('Number of perfect guide matches: ' + str(num_p_matches) + '\n')
                statfile.write('Number of nonperfect guide matches: ' + str(num_np_matches) + '\n')
                statfile.write('Number of undetected guides: ' + str(guides_no_reads) + '\n')
                statfile.write('Percentage of unmapped reads (key not found): ' + str(pct_unmapped) + '\n') #
                statfile.write('Percentage of guides that matched perfectly: ' + str(pct_p_match) + '\n') #
                statfile.write('Percentage of undetected guides: ' + str(pct_no_reads) + '\n') #
                statfile.write('Skew ratio of top 10% to bottom 10%: ' + str(skew_ratio) + '\n') #
                statfile.write('Read coverage: ' + str(coverage))
                statfile.close()
                print(str(in_fastq), 'processed')

    # export files and return dataframes if necessary
    path = Path.cwd()
    if save: 
        outpath = path / out_dir
        Path.mkdir(outpath, exist_ok=True)
        df_ref.to_csv(outpath / out_file, index=False)
        print('count_reads outputed to', str(outpath / out_file))
    print('Count reads completed')
    if return_df:
        return df_ref
        
