"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: }
"""

# imports
import argparse # for parsing arguments from command line

# importing functions
from gene import GeneForCRISPR
from base_editing_guides import find_gRNAs


parser = argparse.ArgumentParser(description='find all guides accessible for base editing')

## all required arguments
parser.add_argument('gene_filepath', type=str,
                    help='gene sequence relative filepath to .fasta file')
parser.add_argument('editing_mode', type=str,
                    help='which base edit is made (ABE CBE) or ')
parser.add_argument('cas_type', type=str,
                    help='which Cas9 and associated PAM is used')
parser.add_argument('protein_filepath', type=str,
                    help='protein sequence relative filepath to .fasta file')

## all optional arguments
parser.add_argument('--guide_length', type=int,
                    help='length of the guide, default is 23')
parser.add_argument('--output_dir', type=str,
                    help='output directory for the guide .csv file')
parser.add_argument('--window', type=str,
                    help='the editing window inclusive entered as , default is 4 to 8')
parser.add_argument('--PAM', type=str,
                    help='a motif which is required for recognition of a guide, entering this overrides the cas_type argument')


# parse arguments from command line
args = parser.parse_args()
print(args)

# create gene object and parses all guides as preprocessing
gene = GeneForCRISPR(filepath=args.gene_filepath)
print('Successfully created gene object from', args.gene_filepath)
gene.parse_exons()
print('Successfully parsed exons')
gene.find_all_guides()
print('Successfully parsed guides')

find_gRNAs(aa_file     = args.protein_filepath,
           gene_object = gene, 
           mode        = args.editing_mode, 
           cas_type    = args.cas_type
          ).to_csv(args.cas_type+args.editing_mode+'_library.csv')
# '230524_FOXA1_BE'+

print('Successfully identified guides')
print('Successfully annotated guides')


# split up find_gRNAs into filter and annotate functions
# add an option to change the output directory of the final .csv, add a prefix to the output filename
# 
