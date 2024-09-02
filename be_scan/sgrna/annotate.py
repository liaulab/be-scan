"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate window and mutation information about the guide}
"""

import pandas as pd
import numpy as np
from pathlib import Path

from _genomic_ import *
from _guideRNA_ import *

def annotate(guides_file, edit_from, edit_to,

    protein_filepath='', window=[4,8], 
    seq_col = 'sgRNA_seq', gene_pos_col='gene_pos', frame_col = 'starting_frame', 
    strand_col = 'sgRNA_strand', window_start_col='windowstart_pos', window_end_col='windowend_pos', 
    output_name="annotated.csv", output_dir='',
    return_df=True, save_df=True,
    ): 
    
    """[Summary]
    Annotates a list of guides in a dataframe with mutational information. 

    Parameters
    ------------
    guides_file: str or path
        The file with the list of guide sequences
    edit_from: str
        The base (ACTG) to be replaced, can be a string of multiple bases
    edit_to: str
        The base (ACTG) to replace with, can be a string of multiple bases

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
       'sgRNA_ID'       : str,    the ID code for the guide
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
    col_names = frame_col, strand_col, gene_pos_col, seq_col, window_start_col, window_end_col

    # READ GUIDES FILE #
    guides_filepath = Path(guides_file)
    df = pd.read_csv(guides_filepath)

    # ASSERTIONS #
    for col in col_names: 
        assert col in df.columns, f"Error {col} not found"

    if len(protein_filepath) > 0: amino_acid_seq = protein_to_AAseq(protein_filepath)
    else: amino_acid_seq = None

    # sgRNA_ID #
    if 'sgRNA_ID' not in df: 
        df.insert(loc=0, column='sgRNA_ID', value=['sgRNA_'+str(i) for i in range(len(df))])

    edit = edit_from, edit_to
    pre = edit_from + 'to' + edit_to
    # coding_seq, THE CODING SEQUENCE IN THE GENOME BEING EDITED #
    df['coding_seq'] = np.where(df[strand_col]=='sense', df[seq_col], 
                                df[seq_col].apply(lambda x: rev_complement(complements, x)) )
    # DELETE ENTRIES WITH THE SAME CODING SEQUENCE #
    dupl_rows = df.duplicated(subset='coding_seq', keep=False)
    df = df[~dupl_rows]

    # win_overlap #
    df[f'{pre}_win_overlap'] = df['coding_seq'].apply(lambda x: annotate_intron_exon(x, window))
    # edit_from+'_count'
    df[edit_from+'_count'] = df['sgRNA_seq'].apply(lambda x: sum([x[window[0]-1:window[1]].count(e) for e in edit_from]))

    # CALCULATE target_CDS, codon_window, residue_window #
    df[f'{pre}_target_CDS'] = df.apply(lambda x: calc_coding_window(x, window, col_names), axis=1) 
    df[[f'{pre}_codon_window', f'{pre}_residue_window']] = df.apply(lambda x: calc_target(x, window, col_names), axis=1, result_type='expand')
    # df[f'{pre}_edit_site'] = df[f'{pre}_editing_window'].apply(lambda x: None if x is None else ((x[0]+x[1])/2)//3)

    # PREDICT POSSIBLE MUTATIONS #
    if len(edit_from) > 1: df[f'{pre}_mutations'] = df.apply(lambda x : annotate_dual_muts(x, edit, amino_acid_seq, col_names, pre, window), axis=1)
    elif len(edit_from) == 1: df[f'{pre}_mutations'] = df.apply(lambda x : annotate_muts(x, edit, amino_acid_seq, col_names, pre, window), axis=1)

    # CALC muttypes LIST AND SINGLE muttype #
    df[f'{pre}_muttypes'] = df.apply(lambda x: determine_mutations(x, col_names, pre), axis=1)
    df[f'{pre}_muttype'] = df.apply(lambda x: categorize_mutations(x, pre), axis=1)

    # DROP UNNECESSARY COLUMNS #
    df = df.drop([f'{pre}_target_CDS', f'{pre}_codon_window', f'{pre}_residue_window'], axis=1)

    print(f'Guides annotated for {edit_from} to {edit_to}.')
    if save_df: 
        Path.mkdir(path / output_dir, exist_ok=True)
        df.to_csv(path / output_dir / output_name, index=False)
    if return_df: return df
