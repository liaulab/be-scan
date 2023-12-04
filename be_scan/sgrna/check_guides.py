"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: filtering }
"""

import pandas as pd
import numpy as np

from be_scan.sgrna._genomic_ import rev_complement, complements
import ahocorasick # https://github.com/WojciechMula/pyahocorasick

def check_guides(guides_file,
                 genome_file, 

                 output_name="output.csv",
                 output_dir='',
                ): 
    """
    Generates a list of guides based on a gene .fasta file,
    and filtering these guides based on PAM and edit available
    in a given window. 

    Parameters
    ------------
    guides_file: str or path
        The file with the gene .fasta sequence
    genome_file: str or path
        The file with the genome .fasta sequence

    output_name : str or path, default 'guides.csv'
        Name of the output .csv guides file
    output_dir : str or path, defailt ''
        Directory path of the output .cs guides file

    Returns: 
    ------------
    merged_df : pandas dataframe
        Contains fwd and rev guides in 'sgRNA_seq', and columns 
        'starting_frame', 'sgRNA_pos', 'exon', 'sgRNA_strand', 
        'gene_strand', 'editing_window', 'gene', 'coding_seq'
        and also contains 'genome_occurrences'
    """
    
    # import guides df
    guides_df = pd.read_csv(guides_file)
    # coding_seq is what needs to be input into algorithm since ref is coding
    if 'coding_seq' not in guides_df.columns: 
        guides_df['coding_seq'] = np.where(guides_df['sgRNA_strand'] == 'sense', 
                                           guides_df['sgRNA_seq'], 
                                           guides_df['sgRNA_seq'].apply(lambda x: rev_complement(complements, 
                                                                                                 x)))
    guides_list = list(guides_df.coding_seq)
    # makes dictionary of guides and how often they come up in the genome
    guides_dict = dict(zip(guides_list, [0]*len(guides_list)))
    
    # create automata and load in guides and make automaton
    automaton = ahocorasick.Automaton()
    for idx, key in enumerate(guides_list):
       automaton.add_word(key, (idx, key))
    automaton.make_automaton()
    # start automaton iteration
    it = automaton.iter("")
    
    # read in genome file
    genome = open(genome_file, 'r')
    # read genome file in line by line
    count = 0
    while True: 
        # get next line
        line = genome.readline().rstrip()
        if line == 'N'*80: 
            continue

        # end of file is reached
        if len(line) == 0 or not line:
            print(count, 'lines processed from', genome_file)
            break
        
        # update df with each iteration through genome line
        it.set(line, False)
        for _, (_, guide) in it: 
            guides_dict[guide] += 1

        count += 1
    
    # convert dict to df
    counts_df = pd.DataFrame(list(guides_dict.items()), 
                             columns=['coding_seq', 'genome_occurrences']) 
    # merge 2 dfs so that genome_occurrences of each guide is logged   
    merged_df = pd.merge(guides_df, counts_df, on='coding_seq', how='inner')

    # output merged_df
    merged_df.to_csv(output_dir+output_name)
    return merged_df
