"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate window and mutation information about the guide}
"""

import pandas as pd
import numpy as np

from be_scan.sgrna._genomic_ import rev_complement, protein_to_AAseq
from be_scan.sgrna._genomic_ import complements
from be_scan.sgrna._guideRNA_ import calc_target, calc_coding_window, calc_editing_window
from be_scan.sgrna._guideRNA_ import annotate_mutations, categorize_mutations
from be_scan.sgrna._gene_ import GeneForCRISPR

def annotate_guides(guides_file, gene_filepath, protein_filepath,
                    edit_from, edit_to,
                    window=(4,8), 

                    seq_col = 'sgRNA_seq', frame_col = 'starting_frame', strand_col = 'sgRNA_strand',
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
        column name of the dataframe with the guide sequence
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

    # read in guides file
    guides_df = pd.read_csv(guides_file)
    if seq_col not in guides_df.columns: 
        print('Error', seq_col, 'not found')
        return
    if frame_col not in guides_df.columns: 
        print('Warning', frame_col, 'not found')
        if strand_col not in guides_df.columns: 
            print('Error', strand_col, 'not found. No information about direction (sense, antisense).')
            return
        else: 
            guides_df = parse_frames(guides_df, gene_filepath)
    # after this point we should have sgRNA_seq, starting frame, sgRNA_strand, and window
    # read in protein_file
    amino_acid_seq = protein_to_AAseq(protein_filepath)

    # coding_seq
    guides_df['coding_seq'] = np.where(guides_df['sgRNA_strand']=='sense', 
                                       guides_df['sgRNA_seq'],
                                       guides_df['sgRNA_seq'].apply(lambda x: rev_complement(complements, x))
                                       )

    # calculate editing_window
    guides_df['editing_window'] = guides_df.apply(lambda x: calc_editing_window(x, window), axis=1)
    # win_overlap
    guides_df['win_overlap'] = np.where(guides_df['coding_seq'].apply(lambda x: x[window[0]-1:window[1]].isupper()), 
                                        "Exon", "Exon/Intron"
                                        )

    # edit_from+'_count'
    guides_df[str(edit_from)+'_count'] = guides_df['sgRNA_seq'].apply(lambda x: x[window[0]-1:window[1]].count(edit_from))
    # calculate target_CDS
    guides_df['target_CDS'] = guides_df.apply(lambda x: calc_coding_window(x, window), axis=1) 

    # calculate target_residue
    guides_df['codon_window'] = guides_df.apply(lambda x: calc_target(x, window, 'DNA'), axis=1)
    # calculate target_residue
    guides_df['residue_window'] = guides_df.apply(lambda x: calc_target(x, window, 'AA'), axis=1)

    # calculate edit_site
    guides_df['edit_site'] = guides_df['editing_window'].apply(lambda x: ((x[0]+x[1])/2)//3)

    # calculate mutations
    guides_df['mutations'] = guides_df.apply(lambda x : annotate_mutations(x, edit_from, edit_to, amino_acid_seq), axis=1)
    # calculate muttype
    guides_df['muttypes'] = guides_df['mutations'].apply(categorize_mutations)
    # calculate muttype
    guides_df['muttype'] = guides_df['muttypes'].apply(lambda x: x[0] if len(x) == 1 else 'Mixed')

    # save df
    if save_df: 
        guides_df.to_csv(output_dir+output_name, index=False)
    if return_df: 
        return guides_df

def parse_frames(guides_df, gene_filepath): 

    # alternate route for finding guide frame since PAM/edit/window info not provided for this function
    gene = GeneForCRISPR(filepath=gene_filepath)
    gene.parse_exons()
    gene_seq = ''.join(gene.exons)

    return guides_df

