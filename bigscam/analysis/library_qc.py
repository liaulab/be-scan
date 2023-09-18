"""
Author: Calvin Xiao Yang Hu
Date: 230916

{Description: primary workflow for QC on library validation after library cloning and running NGS}
"""

from _library_ import count_spacers_fwd
import argparse

parser = argparse.ArgumentParser(description='Analyze sequencing data for sgRNA library distribution')

parser.add_argument('-f', '--fastq', type=str, dest='in_fastq',
                    help='fastq file name')
parser.add_argument('-i', '--input', type=str, dest='in_lib',
                    help='input file name')
parser.add_argument('-o', '--output_save', type=str, dest='out',
                    help='output file name, if no parameters are given, the file will not be saved', default=None)
                    
args = parser.parse_args()
results = count_spacers_fwd(args.in_fastq, args.in_lib, args.out)

# lib_size, num_reads, key_not_found, perfect_matches, np_matches, 
# guides_no_reads, pct_perfmatch, pct_no_reads, undetected_str, skew_ratio
print('Size of library: ' + str(results[0]))
print('Number of reads processed: ' + str(results[1]))
print('Number of reads where key was not found: ' + str(results[2]))
print('Number of perfect guide matches: ' + str(results[3]))
print('Number of nonperfect guide matches: ' + str(results[4]))
print('Number of undetected guides: ' + str(results[5]))
print('Percentage of guides that matched perfectly: ' + str(results[6]))
print('Percentage of undetected guides: ' + str(results[7]))
print('List of undetected gRNAs: ' + results[8])
print('Skew ratio of top 10% to bottom 10%: ' + str(results[9]))

# add parameters for the KEY_REGION_START, KEY_REGION_END, KEY
