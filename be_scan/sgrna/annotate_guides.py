"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate window and mutation information about the guide}
"""

import pandas as pd
import numpy as np
from pathlib import Path

from be_scan.sgrna._genomic_ import *
from be_scan.sgrna._guideRNA_ import *

def annotate_guides(guides_file,
                    edit_from, edit_to,
                    
                    domains={}, 
                    gene_filepath='', protein_filepath='', window=[4,8], 
                    seq_col = 'sgRNA_seq', gene_pos_col='gene_pos',
                    frame_col = 'starting_frame', strand_col = 'sgRNA_strand',
                    output_name="annotated.csv", output_dir='',
                    return_df=True, save_df=True,
                   ): 
    """
    Annotates a list of guides in a dataframe with mutational information. 

    Parameters
    ------------
    guides_file: str or path
        The file with the list of guide sequences
    edit_from: str
        The base (ACTG) to be replaced, can be a string of multiple bases
    edit_to: str
        The base (ACTG) to replace with, can be a string of multiple bases

    domains: dict, default {}
        A dictionary of {'domain_name':(range as a tuple)}
    gene_filepath: str or path, default ''
        The file with the gene .fasta sequence
    protein_filepath: str or path, default ''
        The file with the protein .fasta sequence for double checking the mutations annotated
    window: tuple or list, default = [4,8]
        Editing window, 4th to 8th bases inclusive by default

    seq_col : str, default 'sgRNA_seq'
        column name of the dataframe with the guide sequence, 
        this column must exist
    frame_col : str, default 'starting_frame'
        column name of the dataframe with the frame of the 
        first nucleotide of the guide (0, 1, 2)
    strand_col : str, default 'sgRNA_strand'
        column name of the dataframe with the strand direction (sense, antisense)
    output_name : str or path, default 'annotated.csv'
        Name of the output .csv guides file
    output_dir : str or path, default ''
        Directory path of the output .cs guides file
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save_df : bool, default True
        Whether or not to save the resulting dataframe

    Returns
    ------------
    df : pandas dataframe
    Dataframe contains
       'coding_seq'     : str,    the sense strand sequence of the guide, always fwd
       'editing_window' : tuple,  the gene positions of the editing windows bounds inclusive
       'win_overlap'    : str,    where the window sits (Exon, Exon/Intron, Intron)
       'X_count'        : int,    number of nucleotides in window of A/C/G/T in guide
       'target_CDS'     : str,    the coding nucleotides that are in the editing window
       'codon_window'   : str,    the amino acids that are correspond to the editing window
       'residue_window' : str,    the amino acids that are correspond to the editing window
       'edit_site'      : int,    the amino acid index corresponding to middle nucleotide in editing window
       'mutations'      : str,    a list of mutation (ie F877L, F877P, F877L/F877P)
       'muttypes'       : list,   Missense Nonsense Silent No_C/Exon EssentialSpliceSite Control unique list
       'muttype'        : str,    muttypes condensed down to one type
    """

    path = Path.cwd()
    col_names = frame_col, strand_col, gene_pos_col, seq_col

    # checks editing information is correct, and if len > 1 it is a dual/multi editor
    assert len(edit_from) == len(edit_to)
    for i in range(len(edit_from)): 
        assert edit_from[i] in bases and edit_to[i] in bases 
    edit = edit_from, edit_to

    # read in guides file
    guides_filepath = Path(guides_file)
    guides_df = pd.read_csv(guides_filepath)

    # assertions and warnings
    assert seq_col in guides_df.columns, f"Error {seq_col} not found"
    if frame_col not in guides_df.columns: 
        print('Warning', frame_col, 'not found')
    if strand_col not in guides_df.columns: 
        print('Warning', strand_col, 'not found. No information about direction (sense, antisense).')

    if not all(name in guides_df.columns for name in col_names):
        guides_df = parse_position_frames(guides_df, seq_col, gene_filepath, 
                                              frame_col, strand_col, gene_pos_col)
    # after this point we should have sgRNA_seq, starting frame, sgRNA_strand, and window
    # read in protein_file
    if len(protein_filepath) > 0: 
        amino_acid_seq = protein_to_AAseq(protein_filepath)
    else: 
        amino_acid_seq = None

    for name in col_names: 
        assert name in guides_df.columns 
    edit = edit_from, edit_to
    prefix = edit_from + 'to' + edit_to
    # coding_seq
    guides_df['coding_seq'] = np.where(guides_df['sgRNA_strand']=='sense', 
                                       guides_df['sgRNA_seq'],
                                       guides_df['sgRNA_seq'].apply(lambda x: rev_complement(complements, x))
                                       )
    # delete entries with duplicates between fwd, between rev, and across fwd and rev
    dupl_rows = guides_df.duplicated(subset='sgRNA_seq', keep=False)
    guides_df = guides_df[~dupl_rows]
    # dupl_rows = guides_df.duplicated(subset='coding_seq', keep=False)
    # guides_df = guides_df[~dupl_rows]

    # calculate editing_window
    guides_df[prefix+'_editing_window'] = guides_df.apply(lambda x: calc_editing_window(x, window, col_names), axis=1)
    # win_overlap
    guides_df[prefix+'_win_overlap'] = guides_df['coding_seq'].apply(lambda x: annotate_intron_exon(x, window))

    # edit_from+'_count'
    if len(edit_from) > 1 and len(edit_to) > 1: 
        guides_df[edit_from+'_count'] = guides_df['sgRNA_seq'].apply(lambda x: sum([x[window[0]-1:window[1]].count(e) for e in edit_from]))
    else: 
        guides_df[edit_from+'_count'] = guides_df['sgRNA_seq'].apply(lambda x: x[window[0]-1:window[1]].count(edit_from))

    # calculate target_CDS
    guides_df[prefix+'_target_CDS'] = guides_df.apply(lambda x: calc_coding_window(x, window, col_names), axis=1) 
    # calculate target_residue
    guides_df[prefix+'_codon_window'] = guides_df.apply(lambda x: calc_target(x, window, 'DNA', col_names), axis=1)
    # calculate target_residue
    guides_df[prefix+'_residue_window'] = guides_df.apply(lambda x: calc_target(x, window, 'AA', col_names), axis=1)
    # calculate edit_site
    guides_df[prefix+'_edit_site'] = guides_df[prefix+'_editing_window'].apply(lambda x: None if x is None else ((x[0]+x[1])/2)//3)

    # calculate mutations
    if len(edit_from) > 1 and len(edit_to) > 1: 
        guides_df[prefix+'_mutations'] = guides_df.apply(lambda x : annotate_dual_mutations(x, edit, amino_acid_seq, col_names, prefix), axis=1)
    else:
        guides_df[prefix+'_mutations'] = guides_df.apply(lambda x : annotate_mutations(x, edit, amino_acid_seq, col_names, prefix), axis=1)

    # calculate muttype
    guides_df[prefix+'_muttypes'] = guides_df.apply(lambda x: categorize_mutations(x, col_names, prefix), axis=1)
    # calculate muttype
    guides_df[prefix+'_muttype'] = guides_df[prefix+'_muttypes'].apply(lambda x: None if x is None else x[0] if len(x) == 1 else 'Mixed')

    print('Guides annotated')
    # save df
    if save_df: 
        outpath = path / output_dir
        Path.mkdir(outpath, exist_ok=True)
        guides_df.to_csv(outpath / output_name, index=False)
    if return_df: 
        return guides_df
        