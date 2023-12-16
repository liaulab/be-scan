"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: Generates a dataframe of guides based on taking a gene file and filtering conditions}
"""

import pandas as pd

from be_scan.sgrna._genomic_ import bases, cas_key
from be_scan.sgrna._guideRNA_ import filter_guide, filter_repeats
from be_scan.sgrna._gene_ import GeneForCRISPR
from be_scan.sgrna._genomic_ import process_PAM, rev_complement, complements

def generate_BE_guides(gene_filepath, gene_name, 
                       cas_type, edit_from, edit_to, 
                       PAM=None, window=[4,8], 
                       output_name='guides.csv', output_dir='',
                       return_df=True, save_df=True,
                       ): 
    """
    Generates a list of guides based on a gene .fasta file,
    and filtering these guides based on PAM and edit available
    in a given window. 

    Parameters
    ------------
    gene_filepath: str or path
        The file with the gene .fasta sequence
    gene_name: str
        The name of the gene, can be any string
    cas_type: str
        A type of predetermined Cas (ie Sp, SpG, SpRY, etc)
        This variable is superceded by PAM
    edit_from: char
        The base (ACTG) to be replaced
    edit_to: char
        The base (ACTG) to replace with
    PAM: str, default None
        Optional field to input a custom PAM or a known PAM
        This field supercedes cas_type
    window: tuple or list, default = (4,8)
        Editing window, 4th to 8th bases inclusive by default

    output_name : str or path, default 'guides.csv'
        Name of the output .csv guides file
    output_dir : str or path, default ''
        Directory path of the output .cs guides file
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save_df : bool, default True
        Whether or not to save the resulting dataframe

    Returns
    ------------
    df_no_duplicates : pandas dataframe
    Dataframe contains: 
       'sgRNA_seq'      : str,    the sequence of the guide fwd if on sense strand and rev if on antisense
       'starting_frame' : int,    (0, 1, 2) coding position of the first bp in fwd sgRNA or last bp in rev sgRNA
       'chr_pos'        : int,    the genome position of the first bp in a fwd sgRNA or last bp of a rev sgRNA
       'gene_pos'       : int,    the gene position of the first bp in a fwd sgRNA or last bp of a rev sgRNA
       'exon'           : int,    the exon number according to the input gene_file
       'sgRNA_strand'   : str,    (ie sense or antisense)
       'gene_strand'    : str,    (ie plus or minus)
       'gene'           : str,    name of the gene
    """
    
    # create gene object and parses all guides as preprocessing
    gene = GeneForCRISPR(filepath=gene_filepath)
    print('Create gene object from', gene_filepath)
    gene.parse_exons()
    print('Parsing exons:', len(gene.exons), 'exons found')
    gene.extract_metadata()
    gene.find_all_guides()
    print('Preprocessing sucessful')

    # process cas_type
    if cas_type not in list(cas_key.keys()): 
        raise Exception('Improper cas type input, the options are '+str(list(cas_key.keys())))
    
    # checks editing information is correct
    assert edit_from in bases and edit_to in bases
    edit = edit_from, edit_to
    
    # process PAM, PAM input overrides cas_type
    if PAM is None: 
        PAM = cas_key[cas_type]
    PAM_regex = process_PAM(PAM)

    # process window
    assert window[1] >= window[0] and window[0] >= 0 
    assert window[1] <= len(gene.fwd_guides[0][0])

    # filter for PAM and contains editable base in window
    #    (seq, frame012 of first base, index of first base, exon number)
    fwd_results = [g.copy() for g in gene.fwd_guides if filter_guide(g, PAM_regex, PAM, edit, window)]
    rev_results = [g.copy() for g in gene.rev_guides if filter_guide(g, PAM_regex, PAM, edit, window)]

    # filter out repeating guides in fwd_results rev_results list
    fwd_results = filter_repeats(fwd_results)
    rev_results = filter_repeats(rev_results)

    # adding extra annotations for fwd and rev
    for x in fwd_results: 
        x.append(x[0])
        x.append('sense')
        x.append(gene.strand)
        x.append(gene_name)
    for x in rev_results: 
        x.append(rev_complement(complements, x[0]))
        x.append('antisense')
        x.append(gene.strand)
        x.append(gene_name)

    # set column names for outputing dataframe
    column_names = ['sgRNA_seq', 'starting_frame', 'gene_pos', 'chr_pos', 'exon', 
                    'coding_seq', 'sgRNA_strand', 'gene_strand', 'gene', 
                    ]

    # delete entries with duplicates between fwd, between rev, and across fwd and rev
    df = pd.DataFrame(fwd_results + rev_results, columns=column_names)
    dupl_rows = df.duplicated(subset='sgRNA_seq', keep=False)
    df = df[~dupl_rows]
    dupl_rows = df.duplicated(subset='coding_seq', keep=False)
    df = df[~dupl_rows]
    df_filt = df[~df['sgRNA_seq'].isin(df['coding_seq']) & ~df['coding_seq'].isin(df['sgRNA_seq'])]

    print('Guides generated and duplicates removed')
    # output df
    if save_df: 
        df_filt.to_csv(output_dir + output_name, index=False)
    if return_df: 
        return df_filt
