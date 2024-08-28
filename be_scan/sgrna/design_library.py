"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate window and mutation information about the guide}
"""
import os
from pathlib import Path

from generate_library import generate_library
from reference_check import reference_check
from annotate import annotate

def design_library(gene_filepath, 
                   cas_type, edit_from, edit_to, 

    genome_file='', protein_filepath='', 
    gene_name='', PAM=None, window=[4,8], 

    output_name='annotated_guides.csv', output_dir='', delete=False,
    return_df=True, save_df=True,
    ): 
    
    """[Summary]
    Generates a list of guides based on a gene .fasta file,
    and filtering these guides based on PAM and edit available
    in a given window. 

    Parameters
    ------------
    gene_filepath: str or path
        The file with the gene sequence
    cas_type: str
        A type of predetermined Cas (ie Sp, SpG, SpRY, etc)
        This variable is superceded by PAM
    edit_from: char
        The base (ACTG) to be replaced
    edit_to: char
        The base (ACTG) to replace with

    protein_filepath: str or path, default ''
        The file with the protein .fasta sequence for double checking the mutations annotated
    genome_file: str or path
        The file with the genome sequence. If no path is provided, the reference will not be checked. 
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
       'sgRNA_strand'   : str,    (ie sense or antisense)
       'gene_strand'    : str,    (ie plus or minus)
       'gene'           : str,    name of the gene
       'domain'         : str,    name of the domain according to input ranges, defulats to 'No Domain'
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
       'genome_occurrences' : int, how many times this sequence occurs in the referecnce genome
    """
    
    temp = "temp.csv"

    # GENERATE LIBRARY #
    generate_library_params = {
        'gene_filepath':gene_filepath, 'gene_name':gene_name, 
        'cas_type':cas_type, 'edit_from':edit_from, 'edit_to':edit_to, 
        'PAM':PAM, 'window':window, 'return_df':True, 'save_df':False, 
        }
    guides = generate_library(**generate_library_params)
    guides.to_csv(temp, index=False)

    # ANNOTATE LIBRARY #
    annotate_params = {
        'guides_file':temp, 'gene_filepath':gene_filepath, 
        'protein_filepath':protein_filepath, 'edit_from':edit_from, 'edit_to':edit_to,
        'window':window, 'return_df':True, 'save_df':False
        }
    annotated = annotate(**annotate_params)
    annotated.to_csv(temp, index=False)

    # IF GENOME FILE IS PROVIDED, CHECK GUIDES GAIANST THIS REFERENCE SEQUENCE
    if len(genome_file) > 0: 
        reference_check_params = {
            'guides_file':temp, 'genome_file':genome_file, 
            'delete':delete, 'return_df':True, 'save_df':False
            }
        annotated = reference_check(**reference_check_params)

    os.remove(temp)
    print('Complete! Library generated from', str(gene_filepath))

    if save_df: 
        out_filepath = Path(output_dir)
        annotated.to_csv(out_filepath / output_name, index=False)
    if return_df: 
        return annotated
    
design_library(
    gene_filepath='tests/test_data/sgrna/230408_AR_Input.fasta', 
    cas_type='SpG', 
    edit_from='C', 
    edit_to='T', 
)
