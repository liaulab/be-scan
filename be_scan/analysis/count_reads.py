"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: }
"""

from collections import Counter
import gzip
import warnings

import numpy as np
import pandas as pd

# Count the reads in a FASTQ file and assign them to a reference sgRNA set.
def count_reads(in_fastq, in_ref, KEY_INTERVAL=(10,80),
                KEY='CGAAACACC', KEY_REV='GTTTTAGA', out_counts='counts.csv',
                out_np='np_counts.csv', out_stats='stats.txt', dont_trim_G=False):
    """
    Given a set of sgRNA sequences and a FASTQ file, count the reads in the
    FASTQ, assign the reads to sgRNAs, and export the counts to a csv file `out_counts`. All
    sgRNA sequences not found in the reference file (non-perfect matches) are
    written to a separate csv file `out_np`.

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
    KEY : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_REV : str, default 'GTTTTAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    out_counts : str or path, default 'counts.csv'
        String or path for the output csv file with perfect sgRNA matches.
    out_np : str or path, default 'np_counts.csv'
        String or path for the output csv file with non-perfect sgRNA matches.
    out_stats : str or path, default 'stats.txt'
        String or path for the output txt file with the read counting statistics.
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.
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
    if in_fastq.endswith('.gz'):
        handle = gzip.open(in_fastq, 'rt')
    else:
        handle = open(in_fastq, 'rt')

    # STEP 1B: SET UP VARIABLES FOR SCRIPT
    # make dictionary to hold sgRNA counts - sgRNA_seq, count as k,v
    dict_perfects = {sgRNA:0 for sgRNA in df_ref['sgRNA_seq']}
    list_np = [] # placeholder list for non-perfect matches
    num_reads = 0 # total number of reads processed
    num_perfect_matches = 0 # count of reads with a perfect match to library
    num_np_matches = 0 # count of reads without a perfect match to library
    num_nokey = 0 # count of reads where key was not found
    num_badlength = 0 # count of sgRNA reads that aren't 20bp
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
        if guide in dict_perfects:
            dict_perfects[guide] += 1
            num_perfect_matches += 1
        else:
            num_np_matches += 1
            list_np.append(guide)
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
        statfile.write('Number of reads where length was not 20bp: ' + str(num_badlength) + '\n')
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