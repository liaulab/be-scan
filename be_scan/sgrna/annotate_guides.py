"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate window and mutation information about the guide}
"""

import pandas as pd
import numpy as np
import re

from be_scan.sgrna._genomic_ import rev_complement, protein_to_AAseq, bases
from be_scan.sgrna._genomic_ import complements
from be_scan.sgrna._guideRNA_ import calc_target, calc_coding_window, calc_editing_window
from be_scan.sgrna._guideRNA_ import annotate_mutations, categorize_mutations, parse_position_frames

def annotate_guides(guides_file, gene_filepath, protein_filepath,
                    edit_from, edit_to,
                    window=(4,8), 

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
        The file with the gene .fasta sequence
    gene_filepath: str or path
        The file with the gene .fasta sequence
    edit_from: char
        The base (ACTG) to be replaced
    edit_to: char
        The base (ACTG) to replace with
    protein_filepath: str or path
        The file with the protein .fasta sequence
    window: tuple or list, default = (4,8)
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

    col_names = frame_col, strand_col, gene_pos_col, seq_col

    # checks editing information is correct
    assert edit_from in bases and edit_to in bases
    edit = edit_from, edit_to

    # read in guides file
    guides_df = pd.read_csv(guides_file)

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
    amino_acid_seq = protein_to_AAseq(protein_filepath)

    for name in col_names: 
        assert name in guides_df.columns 
    edit = edit_from, edit_to
    # coding_seq
    guides_df['coding_seq'] = np.where(guides_df['sgRNA_strand']=='sense', 
                                       guides_df['sgRNA_seq'],
                                       guides_df['sgRNA_seq'].apply(lambda x: rev_complement(complements, x))
                                       )

    # calculate editing_window
    guides_df['editing_window'] = guides_df.apply(lambda x: calc_editing_window(x, window, col_names), axis=1)
    # win_overlap
    guides_df['win_overlap'] = np.where(guides_df['coding_seq'].apply(lambda x: x[window[0]-1:window[1]].isupper()), 
                                        "Exon", "Exon/Intron"
                                        )

    # edit_from+'_count'
    guides_df[str(edit_from)+'_count'] = guides_df['sgRNA_seq'].apply(lambda x: x[window[0]-1:window[1]].count(edit_from))
    # calculate target_CDS
    guides_df['target_CDS'] = guides_df.apply(lambda x: calc_coding_window(x, window, col_names), axis=1) 

    # calculate target_residue
    guides_df['codon_window'] = guides_df.apply(lambda x: calc_target(x, window, 'DNA', col_names), axis=1)
    # calculate target_residue
    guides_df['residue_window'] = guides_df.apply(lambda x: calc_target(x, window, 'AA', col_names), axis=1)

    # calculate edit_site
    guides_df['edit_site'] = guides_df['editing_window'].apply(lambda x: None if x is None else ((x[0]+x[1])/2)//3)

    # calculate mutations
    guides_df['mutations'] = guides_df.apply(lambda x : annotate_mutations(x, edit, amino_acid_seq, col_names), axis=1)
    # calculate muttype
    guides_df['muttypes'] = guides_df.apply(lambda x: categorize_mutations(x, col_names), axis=1)
    # calculate muttype
    guides_df['muttype'] = guides_df['muttypes'].apply(lambda x: None if x is None else x[0] if len(x) == 1 else 'Mixed')

    # save df
    if save_df: 
        guides_df.to_csv(output_dir+output_name, index=False)
    if return_df: 
        return guides_df
        