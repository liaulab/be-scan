"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate how often a guide appears in the genome reference file}
"""

from pathlib import Path
import pandas as pd
import numpy as np

from be_scan.sgrna._genomic_ import rev_complement, complements
import ahocorasick # https://github.com/WojciechMula/pyahocorasick

def reference_check(guides_file, genome_file, 

    output_name="filtered.csv", output_dir='', delete=False,
    return_df=True, save_df=True, 
    ): 

    """[Summary]
    Annotates a list of guides with a count of how many times,
    a guide and it's rev_complement appears in the reference genome file.

    Requirement: guides_file must have a column named 'sgRNA_seq' and a column named 'sgRNA_strand'

    Parameters
    ------------
    guides_file: str or path
        The file with the gene .fasta sequence
    genome_file: str or path
        The file with the genome .fasta sequence

    delete : bool, default False
        Whether or not to delete guides with multiple genomic occurrences
    output_name : str or path, default 'guides.csv'
        Name of the output .csv guides file
    output_dir : str or path, default ''
        Directory path of the output .cs guides file
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save_df : bool, default True
        Whether or not to save the resulting dataframe

    Returns: 
    ------------
    merged_df : pandas dataframe
        Contains fwd and rev guides in 'sgRNA_seq', and columns 
        'starting_frame', 'sgRNA_pos', 'exon', 'sgRNA_strand', 
        'gene_strand', 'editing_window', 'gene', 'coding_seq'
        and also contains 'genome_occurrences'
    """
    
    path = Path.cwd()
    guides_filepath = Path(guides_file)
    # import guides df
    guides_df = pd.read_csv(guides_filepath)
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
    genome_filepath = Path(genome_file)
    genome = open(genome_filepath, 'r')
    # read genome file in line by line
    count = 0
    while True: 
        # get next line
        line = genome.readline().rstrip()
        count += 1
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

    # convert dict to df
    counts_df = pd.DataFrame(list(guides_dict.items()), 
                             columns=['coding_seq', 'genome_occurrences']) 
    # merge 2 dfs so that genome_occurrences of each guide is logged   
    merged_df = pd.merge(guides_df, counts_df, on='coding_seq', how='inner')
    # if user wants, delete guides with logged multiple occurrences in reference genome
    if delete: 
        merged_df = merged_df.drop(merged_df[merged_df.genome_occurrences >= 2].index)

    count = (merged_df['genome_occurrences'] > 1).sum()
    print(count, 'guides out of', merged_df.shape[0], 'guides occurred more than once in the reference genome')
    print('Guides checked against reference genome')
    # output merged_df
    if save_df: 
        outpath = path / output_dir
        Path.mkdir(outpath, exist_ok=True)
        merged_df.to_csv(outpath / output_name, index=False)
    if return_df: 
        return merged_df
