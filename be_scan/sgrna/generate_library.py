"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: Generates a dataframe of guides based on taking a gene file and filtering conditions}
"""

from pathlib import Path
import pandas as pd
import warnings

from be_scan.sgrna._genomic_ import *
from be_scan.sgrna._guideRNA_ import filter_guide, filter_sequence
from be_scan.sgrna._gene_ import GeneForCRISPR
# from _genomic_ import *
# from _guideRNA_ import filter_guide, filter_sequence
# from _gene_ import GeneForCRISPR

def generate_library(gene_filepath, 
                     cas_type, edit_from, edit_to, 

    gene_name='', PAM=None, window=[4,8], 
    output_name='guides.csv', output_dir='', return_df=True, save_df=True, 
    exclude_introns=False, exclude_nonediting=False, exclude_duplicates=True, exclude_sequences=['TTTT']
    ): 
    
    """[Summary]
    Generates a list of guides based on a gene .fasta file,
    and filtering these guides based on PAM and edit available
    in a given window. 

    Parameters
    ------------
    gene_filepath: str or path
        The file with the gene .fasta sequence
    cas_type: str
        A type of predetermined Cas (ie Sp, SpG, SpRY, etc)
        This variable is superceded by PAM
    edit_from: char
        The base (ACTG) to be replaced
    edit_to: char
        The base (ACTG) to replace with

    gene_name: str, default ''
        The name of the gene, can be any string
    PAM: str, default None
        Optional field to input a custom PAM or a known PAM
        This field supercedes cas_type
    window: tuple or list, default = [4,8]
        Editing window, 4th to 8th bases inclusive by default

    output_name : str or path, default 'guides.csv'
        Name of the output .csv guides file
    output_dir : str or path, default ''
        Directory path of the output .cs guides file
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save_df : bool, default True
        Whether or not to save the resulting dataframe
    exclude_introns : bool, default True
        Whether or not the editible base needs to be in an intron
    exclude_nonediting : bool, default True
        Whether or not the editible base needs to be in the window
    exclude_duplicates : bool, default True
        Whether or not duplicate guides should be removed from the pool
    exclude_sequences : list of strings, defailt ['TTTT']
        Exclude guides with sequences in this list

    Returns
    ------------
    df_no_duplicates : pandas dataframe
    Dataframe contains: 
       'sgRNA_seq'      : str,    the sequence of the guide 20 bps fwd if on sense strand and rev if on antisense
       'PAM_seq'        : str,    the sequence of the PAM 3 bps fwd if on sense strand and rev if on antisense
       'starting_frame' : int,    (0, 1, 2) coding position of the first bp in fwd sgRNA or last bp in rev sgRNA
       'chr_pos'        : int,    the genome position of the first bp in a fwd sgRNA or last bp of a rev sgRNA
       'gene_pos'       : int,    the gene position of the first bp in a fwd sgRNA or last bp of a rev sgRNA
       'exon'           : int,    the exon number according to the input gene_file
       'windowstart_pos': int,    the gene position of the first bp in the window
       'windowend_pos'  : int,    the gene position of the last bp in the window
       'sgRNA_strand'   : str,    (ie sense or antisense)
       'gene_strand'    : str,    (ie plus or minus)
       'gene'           : str,    name of the gene
    """

    # PREPROCESS cas_type #
    if cas_type not in list(cas_key.keys()): 
        raise Exception('Improper cas type input, the options are '+str(list(cas_key.keys())))
    # PREPROCESS edit_from edit_to, 2 CHARS INDICATES DUAL EDITOR #
    assert len(edit_from) == len(edit_to)
    for i in range(len(edit_from)): 
        assert edit_from[i] in bases and edit_to[i] in bases
        if edit_from[i] == edit_to[i]: 
            warnings.warn(f'You are mutating from {edit_from[i]} to {edit_to[i]}')
    edit = edit_from, edit_to
    # PREPROCESS pam, pam OVERRIDES cas_type #
    if PAM is None: 
        PAM = cas_key[cas_type]
    PAM_regex = process_PAM(PAM)
    # PREPROCESS WINDOW #
    assert window[1] >= window[0] and window[0] >= 0, "Input a valid window ie [4,8]"
    assert window[1] <= 20, "Input a valid window ie [4,8]"
    
    path = Path.cwd()
    # CREATE GENE OBJECT #
    gene = GeneForCRISPR(filepath=gene_filepath)
    print('Create gene object from', gene_filepath)
    gene.parse_exons()
    print('Parsing exons:', len(gene.exons), 'exons found')
    gene.extract_metadata()
    # PARSE ALL GUIDES #
    gene.find_all_guides(window=window)
    print('Preprocessing sucessful!')
    
    # SET COLUMN NAMES FOR OUTPUT #
    column_names = ['sgRNA_seq', 'PAM_seq', 'starting_frame', 'gene_pos', 'chr_pos', 'exon', 
                    'windowstart_pos', 'windowend_pos', 'UTR', 'coding_seq', 'sgRNA_strand', 
                    'gene_strand', 'gene', ]

    # FILTER LIBRARY ACCORDING TO SPECIFICATIONS FOR NONEDITING, INTRONIC, PAM #
    filter_guide_input = {'PAM_regex':PAM_regex, 'edit':edit, 'window':window, 
                          'excl_introns':exclude_introns, 'excl_nonediting':exclude_nonediting}
    fwd_results = [g.copy() for g in gene.fwd_guides if filter_guide(g, **filter_guide_input)]
    rev_results = [g.copy() for g in gene.rev_guides if filter_guide(g, **filter_guide_input)]
    # FILTER OUT UNWANTED SEQUENCES #
    for sequence in exclude_sequences: 
        fwd_results = filter_sequence(fwd_results, sequence)
        rev_results = filter_sequence(rev_results, sequence)

    # ADD EXTRA ANNOTATIONS AND COMBINE #
    for x in fwd_results: 
        x.append(x[0])
        x.append('sense')
    for x in rev_results: 
        x.append(rev_complement(complements, x[0]))
        x.append('antisense')
    results = fwd_results + rev_results
    for x in results: 
        x.append(gene.strand)
        x.append(gene_name)

    # DELETE DUPLICATES BETWEEN FWD, BETWEEN REV, BETWEEN FWD AND REV #
    df = pd.DataFrame(results, columns=column_names)
    if exclude_duplicates: 
        dupl_rows = df.duplicated(subset='sgRNA_seq', keep=False)
        df = df[~dupl_rows]

    print('Guides generated and duplicates removed')
    print(df.shape[0], 'guides were generated')
    # SAVE AND OUTPUT #
    if save_df: 
        Path.mkdir(path / output_dir, exist_ok=True)
        df.to_csv(path / output_dir / output_name, index=False)
    if return_df: 
        return df, gene
