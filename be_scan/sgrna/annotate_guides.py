"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate how often a guide appears in the genome reference file}
"""

import pandas as pd
import numpy as np

from be_scan.sgrna._genomic_ import rev_complement, protein_to_AAseq
from be_scan.sgrna._genomic_ import complements
from be_scan.sgrna._guides_ import calc_target, calc_coding_window, calc_editing_window
from be_scan.sgrna._guides_ import annotate_mutations, categorize_mutations
from be_scan.sgrna._gene_ import GeneForCRISPR

def annotate_guides(guides_file,
                    gene_filepath,
                    edit_from,
                    edit_to,
                    protein_filepath,

                    sgRNA_seq_col = 'sgRNA_seq',
                    starting_frame_col = 'starting_frame',
                    sgRNA_strand_col = 'sgRNA_strand',
                    window=(4,8), 
                    output_name="annotated.csv",
                    output_dir='',
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
       'muttypes'       : list,   Missense Nonsense Silent No_C/Exon EssentialSpliceSite Control
       'muttype'        : str,    muttypes condensed down to unique types
       'muttype'        : str,    muttypes condensed down to one type
    """

    # read in guides file
    guides_df = pd.read_csv(guides_file)
    if sgRNA_seq_col not in guides_df.columns: 
        print('Error', sgRNA_seq_col, 'not found')
        return
    if starting_frame_col not in guides_df.columns: 
        print('Warning', starting_frame_col, 'not found')
        if sgRNA_strand_col not in guides_df.columns: 
            print('Error', sgRNA_strand_col, 'not found. No information about direction (sense, antisense).')
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
    guides_df.to_csv(output_dir+output_name)
    return guides_df

def parse_frames(guides_df, gene_filepath): 

    # alternate route for finding guide frame since PAM/edit/window info not provided for this function
    gene = GeneForCRISPR(filepath=gene_filepath)
    gene.parse_exons()
    gene_seq = ''.join(gene.exons)

    return guides_df

# this is the main function for taking in lists of guides, 
# then annotating all their predicted edits
# def annotate_BE_guides(protein_filepath, fwd_guides, rev_guides, edit_from, edit_to, window=[4,8]): 
#     ### ADD DOCS
#     # Parameters
#     #    protein_filepath: filepath to an amino acid sequence corresponding to gene file
#     #    fwd_guides, rev_guides: generated from identify_guides
#     #    edit_from: the base (ACTG) to be replaced
#     #    edit_to: the base (ACTG) to replace with
#     #    window: editing window, 4th to 8th bases inclusive by default
#     # outputs: a dataframe
    
#     ### target_codons: list of codons that we want to make with our base edit

#     # codon indices, predicted edits made
#     amino_acid_seq = protein_to_AAseq(protein_filepath)
#     num_aa = 3 # this is the num of amino acids we look ahead in our frame

#     for g in fwd_guides: 
#         # mutates all residues according to the mode, every combination of residue mutations
#         original = g[0][:12] # a string
#         guide_window = g[0][window[0]-1:window[1]] # a string
#         mutateds = [g[0][:window[0]-1] + m + g[0][window[1]:12] for m in make_mutations(guide_window, edit_from, edit_to)] # list of strings 

#         # compares the residues to find which amino acids were altered and catalogs them
#         edits, edit_inds = [], [] # lists of lists, of all edits for all possible mutations
#         start = (-1*g[1])+3
#         orig = original[start:start+(num_aa*3)]
#         for m in mutateds: 
#             # for each possible mutation, come up with the list of amino acid changes
#             edit, edit_ind = find_aa_edits_fwd(m, g, start, orig, num_aa, amino_acid_seq)
#             edits.append(edit)
#             edit_inds.append(edit_ind)
#         # append all information to dataframe
#         g.extend([edits, edit_inds, 'fwd'])
#         assert(len(g)) == 7
        
#     for g in rev_guides: 
#         # mutates all residues according to the mode, every combination of residue mutations
#         original = g[0][:12] # a string
#         guide_window = g[0][window[0]-1:window[1]] # a string
#         mutateds = [g[0][:window[0]-1] + m + g[0][window[1]:12] for m in make_mutations(guide_window, edit_from, edit_to)] # a list of strings

#         # compares the residues to find which amino acids were altered and catalogs them
#         edits, edit_inds = [], [] # lists of lists, of all edits for all possible mutations
#         start = g[1]+1
#         orig = rev_complement(complements, original[start:start+(num_aa*3)])
#         for m in mutateds: 
#             # for each possible mutation, come up with the list of amino acid changes
#             edit, edit_ind = find_aa_edits_rev(m, g, start, orig, num_aa, amino_acid_seq)
#             edits.append(edit)
#             edit_inds.append(edit_ind)
#         # append all information to dataframe
#         g.extend([edits, edit_inds, 'rev'])
#         assert(len(g)) == 7
        
# #     print(pd.DataFrame(fwd_guides + rev_guides))
#     return pd.DataFrame(fwd_guides + rev_guides)
