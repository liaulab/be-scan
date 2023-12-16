"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate window and mutation information about the guide}
"""
import os
from pathlib import Path

from be_scan.sgrna.generate_guides import generate_BE_guides
from be_scan.sgrna.check_guides import check_guides
from be_scan.sgrna.annotate_guides import annotate_guides

def guides(gene_filepath, gene_name, genome_file, protein_filepath, 
           cas_type, edit_from, edit_to, 

           PAM=None, window=(4,8), 
           output_name='annotated_guides.csv', output_dir='', delete=False,
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
    genome_file: str or path
        The file with the genome .fasta sequence
    protein_filepath: str or path
        The file with the protein .fasta sequence
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
    temp = "temp.csv"
    gene_filepath, genome_file, protein_filepath = Path(gene_filepath), Path(genome_file), Path(protein_filepath)
    
    guides = generate_BE_guides(gene_filepath=gene_filepath,
                                gene_name=gene_name, 
                                cas_type=cas_type, 
                                edit_from=edit_from, edit_to=edit_to, 
                                PAM=PAM, window=window, 
                                return_df=True, save_df=False
                                )
    guides.to_csv(temp, index=False)

    filtered = check_guides(temp, 
                            genome_file=genome_file, 
                            delete=delete, 
                            return_df=True, save_df=False)
    filtered.to_csv(temp, index=False)

    annotated = annotate_guides(temp, 
                                gene_filepath=gene_filepath, 
                                protein_filepath=protein_filepath,
                                edit_from=edit_from, edit_to=edit_to,
                                window=window, 
                                return_df=True, save_df=False,
                                )

    os.remove(temp)
    print('Complete! Library generated from', str(gene_filepath))

    if save_df: 
        out_filepath = Path(output_dir)
        annotated.to_csv(out_filepath / output_name, index=False)
    if return_df: 
        return annotated
    