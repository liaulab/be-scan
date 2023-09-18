"""
Author: Kevin C Ngan, Calvin Xiao Yang Hu
Date: 230916

{Description: }
"""

from Bio import SeqIO
import csv
from collections import OrderedDict
import numpy as np
import sys
import gzip
from collections import Counter

def count_spacers_fwd(in_fastq, in_lib, out): 
    
    """
    Creates dictionary with guide counts from fastq_file, writes to output_file
    in_lib: csv file with guide sequences (default: library_sequences.csv)
    in_fastq: forward read fastq file (default: NGS.fastq)
    out_counts: prefix for csv output for guide dictionary, csv output for imperfect matches, txt output for statistics
    FORWARD: Use end of hU6 promoter to identify sgRNA sequence
    """

    # STEP 1: SET THE KEY REGION START/END AND KEY SEQUENCE TO FIND

    KEY_REGION_START = 5 # start index of key region
    KEY_REGION_END = 40 # end index of key region
    KEY = "CGAAACACCG" # identifies seq before guide to determine position

    # STEP 2: RESET VARIABLES, OPEN INPUT FILES FOR PROCESSING

    # set all variables to 0
    num_reads = 0 # total number of reads processed
    perfect_matches = 0 # count of reads with a perfect match to library
    np_matches = 0 # count of reads without a perfect match to library
    key_not_found = 0 # count of reads where key was not found

    # open library sequences and start dictionary, set counts to 0
    # dict_perfects: guide sequence as key, guide count as entry
    with open(in_lib) as infile:
        reader = csv.reader(infile)
        dict_perfects = {rows[0]:0 for rows in reader}
            # build dict with each row in reader as key with value = 0
    try:
        handle = gzip.open(in_fastq, "rt")
    except:
        print("could not find fastq file")
        return
    lib_size = len(dict_perfects)

    # STEP 3: PROCESS FASTQ FILE READS AND ADD COUNTS TO DICT

    # process reads in fastq file
    np_list = [] # generate placeholder list for non-perfect matches
    readiter = SeqIO.parse(handle, "fastq")
    for record in readiter: # contains the seq and Qscore etc.
        num_reads += 1
        read_sequence = str.upper(str(record.seq))
        key_region = read_sequence[KEY_REGION_START:KEY_REGION_END]
        key_index = key_region.find(KEY)
        if key_index >= 0:
            start_index = key_index + KEY_REGION_START + len(KEY)
            guide = read_sequence[start_index:(start_index + 20)]
            if guide in dict_perfects:
                dict_perfects[guide] += 1
                perfect_matches += 1
            else:
                np_matches += 1
                np_list.append(guide)
        else:
            key_not_found += 1
    # create dictionary for non-perfect guide matches
    dict_np = Counter(np_list)

    # STEP 4: CALCULATE STATS

    # percentage of guides that matched perfectly
    pct_perfmatch = round(perfect_matches/float(perfect_matches + np_matches) * 100, 1)

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
        
    undetected_list = []
    for key, item in dict_perfects.items():
        if item == 0:
            undetected_list.append(key)
    
    undetected_str = ''
    for gRNA in undetected_list:
        temp = gRNA + ", "
        undetected_str += temp
    n = len(undetected_str)
    undetected_str = undetected_str[0:n-2]

    handle.close()
    results = (lib_size, num_reads, key_not_found, perfect_matches, np_matches, guides_no_reads, pct_perfmatch, pct_no_reads, undetected_str, skew_ratio)
    
    # STEP 5: save results if an out prefix is provided
    if out is not None: 
        save_results(results, dict_perfects, dict_np, out)
    
    return results


def save_results(results, dict_perfects, dict_np, out): 

    # STEP 0: SET OUTPUT FILE NAMES

    out_counts = out + "_library_count.csv"
    out_np = out + "_imperfect_counts.csv"
    out_stats = out + "_library_statistics.txt"

    # STEP 1: GENERATE SUMMARY STATS INTO STAT OUTPUT FILE

    # write analysis statistics to stat_file
    with open(out_stats, 'w') as statfile:
        statfile.write('Number of reads processed: ' + str(results[1]) + '\n')
        statfile.write('Number of reads where key was not found: ' + str(results[2]) + '\n')
        statfile.write('Number of perfect guide matches: ' + str(results[3]) + '\n')
        statfile.write('Number of nonperfect guide matches: ' + str(results[4]) + '\n')
        statfile.write('Number of undetected guides: ' + str(results[5]) + '\n')
        statfile.write('Percentage of guides that matched perfectly: ' + str(results[6]) + '\n')
        statfile.write('Percentage of undetected guides: ' + str(results[7]) + '\n')
        statfile.write('List of undetected gRNAs: ' + results[8] + '\n')
        statfile.write('Skew ratio of top 10% to bottom 10%: ' + str(results[9]))
        statfile.close()
    
    # STEP 2: SORT DICTIONARIES AND GENERATE OUTPUT FILES

    # create ordered dict with guides (A to Z) and respective counts
    # sort perfect matches dictionary first and output to library_count.csv
    perfdict_sorted = OrderedDict(sorted(dict_perfects.items(), key=lambda t: t[0]))
    with open(out_counts, 'w', newline='') as perfectmatches_csv:
        mywriter = csv.writer(perfectmatches_csv, delimiter=',')
        for guide in perfdict_sorted:
            count = perfdict_sorted[guide]
            mywriter.writerow([guide,count])

    # now sort non-perfect matches and create output csv
    nonperfdict_sorted = OrderedDict(sorted(dict_np.items(), key=lambda t: t[0]))
    with open(out_np, 'w', newline='') as nonperfectmatches_csv:
        mywriter = csv.writer(nonperfectmatches_csv, delimiter=',')
        for guide in nonperfdict_sorted:
            count = nonperfdict_sorted[guide]
            mywriter.writerow([guide,count])
