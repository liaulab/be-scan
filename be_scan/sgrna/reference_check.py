"""
Author: Calvin XiaoYang Hu
Date: 231204

{Description: Annotate how often a guide appears in the genome reference file}
"""

from pathlib import Path
import pandas as pd
import numpy as np

from be_scan.sgrna._genomic_ import rev_complement, complements
# from _genomic_ import rev_complement, complements
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

    output_name : str or path, default 'guides.csv'
        Name of the output .csv guides file
    output_dir : str or path, default ''
        Directory path of the output .cs guides file
    delete : bool, default False
        Whether or not to delete guides with multiple genomic occurrences
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save_df : bool, default True
        Whether or not to save the resulting dataframe

    Returns
    ------------
    merged_df : pandas dataframe
        Contains fwd and rev guides in 'sgRNA_seq', 'sgRNA_strand', 'coding_seq'
        and also contains 'ref_occurrences'
    """
    path = Path.cwd()
    
    # READ GUIDES FILE #
    guides_filepath = Path(guides_file)
    df = pd.read_csv(guides_filepath)
    # coding_seq IS THE INPUT INTO ALGORITHM SINCE REFERENCE IS CODING #
    if 'coding_seq' not in df.columns: 
        df['coding_seq'] = np.where(df['sgRNA_strand'] == 'sense', df['sgRNA_seq'], 
                                    df['sgRNA_seq'].apply(lambda x: rev_complement(complements, x)))
    guides_list = list(df.coding_seq)
    # DICT OF {GUIDE:GUIDE COUNT} #
    guides_dict = dict(zip(guides_list, [0]*len(guides_list)))
    
    # CREATE AUTOMATA, LOAD IN GUIDES #
    automaton = ahocorasick.Automaton()
    for idx, key in enumerate(guides_list):
       automaton.add_word(key, (idx, key))
    automaton.make_automaton()
    # START AUTOMATON ITERATION #
    it = automaton.iter("")

    # READ IN REFERENCE FILE LINE BY LINE #
    genome_filepath = Path(genome_file)
    genome = open(genome_filepath, 'r')
    # read genome file in line by line
    count = 0
    while True: 
        # GET NEXT LINE #
        line = genome.readline().rstrip()
        count += 1
        if line == 'N'*80: 
            continue

        # END OF FILE #
        if len(line) == 0 or not line:
            print(count, 'lines processed from', genome_file)
            break
        
        # UPDATE DF WITH EACH LOOP #
        it.set(line, False)
        for _, (_, guide) in it: 
            guides_dict[guide] += 1

    # CONVERT DICT TO DF #
    counts_df = pd.DataFrame(list(guides_dict.items()), columns=['coding_seq', 'ref_occurrences'])
    # MERGE WITH PREVIOUS DF #
    merged_df = pd.merge(df, counts_df, on='coding_seq', how='inner')
    # OPTIONAL TO DELETE GUIDES WITH MULTIPLE OCCURRENCES #
    if delete: 
        merged_df = merged_df.drop(merged_df[merged_df.ref_occurrences > 1].index)

    count = (merged_df['ref_occurrences'] > 1).sum()
    print('Guides checked against reference genome')
    print(count, 'guides out of', merged_df.shape[0], 'occurred more than once in the reference genome')
    # OUTPUT merged_df
    if save_df: 
        Path.mkdir(path / output_dir, exist_ok=True)
        merged_df.to_csv(path / output_dir / output_name, index=False)
    if return_df: 
        return merged_df
