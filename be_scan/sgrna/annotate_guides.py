"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate how often a guide appears in the genome reference file}
"""

import pandas as pd
import numpy as np

from be_scan.sgrna._gene_ import GeneForCRISPR

def annotate_guides(guides_file,
                    target,

                    output_name="output.csv",
                    output_dir='',
                   ): 
    """
    Annotates a list of guides with 
       'X_count'        : int,    number of nucleotides in window of A/C/G/T
       'target_CDS'     : tuple,  the nucleotides that are in the editing window
       'target_residue' : tuple,  the amino acids that are in the editing window
       'mutations'      : str,    a list of mutation (ie F877L)
       'muttype'        : str,    Missense Nonsense Silent No_C/Exon EssentialSpliceSite Control
       'edit_site'      : int,    the amino acid index corresponding to first nucleotide in editing window
    And outputs an annotated dataframe

    Parameters
    ------------
    guides_file: str or path
        The file with the gene .fasta sequence
    gene_filepath: str or path
        The file with the gene .fasta sequence
    target: char
        (ACGT) which is the target of editing
    """

    # read in guides file
    guides_df = pd.read_csv(guides_file)

    # calculate X_count
    guides_df['X_count']

    # calculate target_CDS
    guides_df['target_CDS']

    # calculate target_residue
    guides_df['target_residue']

    # calculate mutations
    guides_df['mutations']

    # calculate muttype
    guides_df['muttype']

    # calculate edit_site
    guides_df['X_count']

