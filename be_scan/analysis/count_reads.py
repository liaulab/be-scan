"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Count the reads in a FASTQ file and assign them to a reference sgRNA set}
"""

from collections import Counter
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def count_reads(sample_sheet, annotated_lib, 
    
    KEY_INTERVAL=(10,80), KEY='CGAAACACC', KEY_REV='GTTTGAGA', dont_trim_G=False,
    in_dir='', out_dir='', out_file='library_counts.csv', 
    sgRNA_seq_col = 'sgRNA_seq', 
    return_df=True, save_files=True, plot_out_type='pdf', 
    ): 
    
    """[Summary]
    Given a set of sgRNA sequences and a FASTQ file, count the reads in the
    FASTQ, assign the reads to sgRNAs, and export the counts to a csv file `out_counts`. All
    sgRNA sequences not found in the reference file (non-perfect matches) are
    written to a separate csv file `out_nc`. 
    count_reads works by reading through a .fastq file and searching for the guide in between
    the subsequences for KEY and KEY_REV within the KEY_INTERVAL. These counts are then recorded.

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
    annotated_lib : str or path
        String or path to the reference file. annotated_lib must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.

    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_REV : str, default 'GTTTGAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.

    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'counts_library.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save_files : bool, default True
        Whether or not to save individual counts, noncounts, and stats files
    plot_out_type : str, optional, defaults to 'pdf'
        file type of figure output
    """
    path = Path.cwd()
    inpath = Path(in_dir)
    outpath = Path(out_dir)

    # LOAD IN SAMPLE SHEET #
    sample_filepath = inpath / sample_sheet
    sample_df = pd.read_csv(sample_filepath)
    for colname in ['fastq_file', 'counts_file', 'noncounts_file', 'stats_file', 'condition']: 
        if colname not in sample_df.columns.tolist():
            raise Exception(f"{annotated_lib} is missing column: {colname}")
    samples = [list(a) for a in zip(sample_df.fastq_file, sample_df.counts_file, sample_df.noncounts_file, 
                                    sample_df.stats_file, sample_df.condition) ]

    # STEP 1: OPEN INPUT FILES FOR PROCESSING, CHECK FOR REQUIRED FORMATTING
    # LOOK FOR sgRNA_seq COLUMN #
    library_filepath = inpath / annotated_lib
    df_ref = pd.read_csv(library_filepath, header=0)
    if sgRNA_seq_col not in df_ref.columns.tolist(): 
        raise Exception(f'{annotated_lib} is missing column: {sgRNA_seq_col}')
    df_ref[sgRNA_seq_col] = df_ref[sgRNA_seq_col].str.upper()

    for fastq, counts, nc, stats, cond in samples: 
        # FASTQ FILE OF READS AND PATH TO ALL OUTPUT FILES #
        in_fastq = inpath / fastq
        out_counts, out_nc, out_stats = outpath / counts, outpath / nc, outpath / stats

        # STEP 1B: SET UP VARIABLES FOR SCRIPT
        dict_p = {sgRNA:0 for sgRNA in df_ref[sgRNA_seq_col]} # dictionary (sgRNA_seq:sgRNA_counts}
        list_np = [] # list for non-perfect matches
        # reads count of: total, perfect match, non perfect match, no key found, not 20bps
        num_reads, num_p_matches, num_np_matches, num_nokey, num_badlength = 0, 0, 0, 0, 0
        KEY_START, KEY_END = KEY_INTERVAL[0], KEY_INTERVAL[1] # set the key interval

        # TRY OPENING FASTQ #
        handle = gzip.open(in_fastq, 'rt') if str(in_fastq).endswith('.gz') else open(in_fastq, 'rt')
        # STEP 2: PROCESS FASTQ FILE READS AND ADD COUNTS TO DICT
        while True: # contains the seq and Qscore etc
            read = handle.readline()
            if not read: break # end of file
            elif read.startswith("@"): # if line is a read header, next line is the sequence
                read = handle.readline()
            else: continue # next loop
            
            num_reads += 1
            read_sequence = str.upper(str(read))
            key_region = read_sequence[KEY_START:KEY_END]
            key_index = key_region.find(KEY)
            key_rev_index = key_region.rfind(KEY_REV)
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

        # STEP 3: SORT PERF MATCH DICTIONARIES AND GENERATE OUTPUT FILES
        df_perfects = pd.DataFrame(data=dict_p.items(), columns=[sgRNA_seq_col, cond])
        if save_files: 
            # OUTPUT PERFECT COUNTS #
            df_perfects.sort_values(by=cond, ascending=False, inplace=True)
            df_perfects.to_csv(out_counts, index=False)

            # df_ref[cond] column histogram #
            weights = np.ones_like(df_perfects[cond]) / len(df_perfects[cond])
            plt.hist(df_perfects[cond], weights=weights, bins=(len(df_perfects[cond])//5)+1)
            out = stats.split('.')[0] + '_histogram.' + plot_out_type
            plt.title(f"Distributions of sgRNA in {fastq}")
            plt.xlabel('Count of sgRNA')
            plt.ylabel('Proportion of sgRNA')
            plt.savefig(outpath / out, format=plot_out_type)
            plt.clf()
            plt.close()

        # ADD MATCHING COUNTS TO LIBRARY DF #
        df_ref = pd.merge(df_ref, df_perfects, on=sgRNA_seq_col, how='outer')
        df_ref[cond] = df_ref[cond].fillna(0)

        # SORT NON-PERFECT MATCHES BY FREQUENCY AND OUTPUT #
        dict_np = Counter(list_np) # use Counter to tally up np matches
        nc_name = nc.split("/")[-1]
        df_ncmatches = pd.DataFrame(data=dict_np.items(), columns=[sgRNA_seq_col, nc_name])
        if save_files:
            df_ncmatches.sort_values(by=nc_name, ascending=False, inplace=True)
            df_ncmatches.to_csv(out_nc, index=False)

        # STEP 4: CALCULATE STATS AND GENERATE STAT OUTPUT FILE
        if save_files:
            # CALCULATE READ COVERAGE #
            num_guides = df_ref[sgRNA_seq_col].shape[0]
            coverage = round(num_reads / num_guides, 1)
            # PERCENTAGE OF GUIDES THAT MATCH PERFECTLY #
            pct_p_match = round(num_p_matches/float(num_p_matches + num_np_matches) * 100, 1)
            # PERCENTAGE OF UNDETECTED GUIDES #
            vals_p = np.fromiter(dict_p.values(), dtype=int)
            guides_no_reads = np.count_nonzero(vals_p==0)
            pct_no_reads = round(guides_no_reads/float(len(dict_p.values())) * 100, 1)
            # SKEW RATIO TOP 10% TO BOTTOM 10% #
            top_10 = np.percentile(list(dict_p.values()), 90)
            bottom_10 = np.percentile(list(dict_p.values()), 10)
            if top_10 != 0 and bottom_10 != 0: skew_ratio = top_10/bottom_10
            else: skew_ratio = 'Not enough perfect matches to determine skew ratio'
            # CALCULATE NUMBER OF UNMAPPED READS #
            pct_unmapped = round((num_nokey / num_reads) * 100, 2)

            # WRITE OUTPUT TO STATS FILE #
            with open(out_stats, 'w') as statfile:
                stats_list = [
                    f'Number of reads processed: {num_reads}', 
                    f'Number of reads where key was not found: {num_nokey}', 
                    f'Number of reads where length was not 20bp: {num_badlength}', 
                    f'Number of perfect guide matches: {num_p_matches}', 
                    f'Number of nonperfect guide matches: {num_np_matches}', 
                    f'Number of undetected guides: {guides_no_reads}', 
                    f'Percentage of unmapped reads (key not found): {pct_unmapped}', 
                    f'Percentage of guides that matched perfectly: {pct_p_match}', 
                    f'Percentage of undetected guides: {pct_no_reads}', 
                    f'Skew ratio of top 10% to bottom 10%: {skew_ratio}', 
                    f'Read coverage: {coverage}', 
                ]
                statfile.write('\n'.join(stats_list))
                print(f'{str(in_fastq)} processed')

    # SAVE DF AND RETURN #
    if save_files: 
        Path.mkdir(path / out_dir, exist_ok=True)
        df_ref.to_csv(path / out_dir / out_file, index=False)
        print('count_reads output to', str(outpath / out_file))
    print('Count reads completed')
    if return_df:
        return df_ref
