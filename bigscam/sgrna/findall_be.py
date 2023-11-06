"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: primary workflow for generating and annotating a set of base editing guides}
"""

# imports
import argparse # for parsing arguments from command line

# importing functions
from ._gene_ import GeneForCRISPR
from ._BE_guides_ import identify_BE_guides, annotate_BE_guides

def add_parser_args(parser):

    ## all required arguments
    parser.add_argument('-g', '--gene_filepath', type=str,
                        help='gene sequence relative filepath to .fasta file')
    parser.add_argument('-e1', '--edit_from', type=str,
                        help='base to be edited')
    parser.add_argument('-e2', '--edit_to', type=str,
                        help='base to edit to')
    parser.add_argument('-c', '--cas_type', type=str,
                        help='which Cas9 (and associated PAM) is used')
    parser.add_argument('-p', '--protein_filepath', type=str,
                        help='protein sequence relative filepath to .fasta file')

    ## all optional arguments
    parser.add_argument('--guide_length', type=int,
                        help='length of the guide, default is 23')
    parser.add_argument('--output_dir', type=str,
                        help='output directory for the guides output file')
    parser.add_argument('--output_prefix', type=str,
                        help='output prefix for the guides output file')
    parser.add_argument('--output_type', type=str,
                        help='output file type for the guides output file (default is .csv)')
    parser.add_argument('--window', nargs="+", type=int, 
                        help='the editing window inclusive entered as 2 integers, default is 4 8')
    parser.add_argument('--PAM', type=str,
                        help='a motif which is required for recognition of a guide, entering this overrides the cas_type argument')

    return parser

def main(args):

    # create gene object and parses all guides as preprocessing
    gene = GeneForCRISPR(filepath=args.gene_filepath)
    print('Create gene object from', args.gene_filepath)
    gene.parse_exons()
    print('Parsing exons:', len(gene.exons), 'exons found')
    gene.find_all_guides()
    print('Preprocessing sucessful')

    # sort guides based on PAM/cas type and editing mode first to identify all possible guides
    if args.PAM is None and args.window is None:
        fwd_res, rev_res, mode = identify_BE_guides(gene, args.cas_type, args.edit_from, args.edit_to)
    elif args.PAM is None and args.window is not None:
        fwd_res, rev_res, mode = identify_BE_guides(gene, args.cas_type, args.edit_from, args.edit_to, window=args.window)
    elif args.PAM is not None and args.window is None:
        fwd_res, rev_res, mode = identify_BE_guides(gene, args.cas_type, args.edit_from, args.edit_to, PAM=args.PAM)
    else: 
        fwd_res, rev_res, mode = identify_BE_guides(gene, args.cas_type, args.edit_from, args.edit_to, PAM=args.PAM, window=args.window)
    print('Identifying guides:', len(fwd_res)+len(rev_res), 'guides found')
    # using the preprocessing data to annotate which amino acid residues are mutated
    df = annotate_BE_guides(args.protein_filepath, fwd_res, rev_res, *mode)
    print('Annotating guides successful')

    # correctly format the output path and filetype of the dataframe and save
    if args.output_type is None: 
        output_type = 'csv'
        output_type = output_type.lower()
    else: 
        output_type = args.output_type.lower()
    assert output_type in ['csv', 'tsv', 'dat', 'xls', 'xlsx', 'json']
    if args.output_dir is None: 
        args.output_dir = ''
    if args.output_prefix is None: 
        args.output_prefix = ''
    out_path = args.output_dir+args.output_prefix+'_'+args.cas_type+'_'+args.edit_from+'to'+args.edit_to+'_library.'+output_type
    df.to_csv(out_path)
    print('Successfully saved 0wO')

    # think about how we can reformat this data to account for every combination of edits that can be made with a single guide
