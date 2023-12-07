"""
Author: Calvin XiaoYang Hu
Date: 231102

{Description: helper functions for processing guides and library data}
"""

from itertools import product
from be_scan.sgrna._genomic_ import complements, rev_complement, DNA_to_AA

# evaluates if guide has PAM and has a residue in window
# returns TRUE or FALSE
def filter_guide(g, PAM_regex, PAM, edit, window): 
    """
    Evaluates if a guide has a PAM and target residue within its window. 

    Parameters
    ------------
    g : str, guide sequence
    PAM_regex : regex str, regez string constructued from PAM
    PAM : str, PAM sequence
    edit : tuple, (edit_from edit_to)
    window : tuple, inclusive range for editing

    Returns
    ------------
    bool True or False
    """
    seq = g[0]
    return PAM_regex.match(seq[-len(PAM):]) and edit[0] in seq[window[0]-1:window[1]]

def filter_repeats(results): 
    """
    Delete duplicate guides
    Duplicate guides are difficult to deconvolute

    Parameters
    ------------
    results : list of lists [sgRNA_seq, starting_frame, gene_pos, chr_pos, exon]

    Returns
    ------------
    results : list of lists, shortened input results
    """
    seqs = [g[0] for g in results]
    dup_seqs = set([x for x in seqs if seqs.count(x) > 1])
    results = [g.copy() for g in results if g[0] not in dup_seqs]
    return results

def annotate_mutations(row, edit, amino_acid_seq, col_names): 
    """
    Come up with list of annotations (ie F877L, F877P, F877L/F877P)
    for each guide. 

    Parameters
    ------------
    row : row of input df
    edit : tuple, (edit_from edit_to)
    amino_acid_seq : dict, amino acid with {index : amino acid}
    col_names : tuple, (starting_frame, sgRNA_strand, gene_pos, seq_col)

    Returns
    ------------
    mutation_details : list, a list of mutations
    """
    # extract relevant data from dataframe
    frame, dir, pos = row[col_names[0]], row[col_names[1]], row[col_names[2]]
    dna_window, dna, aa = row['target_CDS'], row['codon_window'], row['residue_window']

    # starting index of amino acid is different for sense vs anti
    if dir == 'sense': 
        start = int((pos+(-1*frame))/3)+2 # indexed at 1, for index at 0 +3 instead of +6
    else: 
        start = int((pos+(-1*frame)+1)/3)-2

    # add to a list of mutations for each row
    mutation_details = []
    for m in mutation_combos(dna_window, edit, dir): 
        if m == dna_window: # if the strand is unmutated
            continue
        # construct new dna and new aa sequence from one potential mutation
        new_dna = dna.replace(dna_window, m)
        new_aa = DNA_to_AA(new_dna, upper=False)
        # write out which mutations are changed
        mutations = format_mutation(aa, new_aa, start, amino_acid_seq)
        mutation_details.append(mutations)
    return mutation_details

def mutation_combos(guide_window, edit, dir):
    """
    Generates a list of all possible mutations of the guide_window
    
    Parameters
    ------------
    guide_window : str, part of a guide that has length divisible by 3
    edit : tuple, (edit_from edit_to)
    dir : str, (sense or antisense)

    Returns
    ------------
    mutated : list, each item is a mutated guide_window
    """
    if dir == 'antisense': 
        edit = (complements[edit[0]], complements[edit[1]])
    mutated = []
    # Convert input string into a list so we can easily substitute letters
    seq = list(guide_window)
    # Find indices of key letters in seq
    indices = [ i for i, c in enumerate(seq) if c in edit[0]]

    # Generate key letter combinations & place them into the list
    for t in product(edit[0]+edit[1], repeat=len(indices)):
        for i, c in zip(indices, t):
            seq[i] = c
        mutated.append(''.join(seq))
    return mutated

def format_mutation(aa, new_aa, start, amino_acid_seq): 
    """
    Translates and formats the mutation with unedited-position-edited
    
    Parameters
    ------------
    aa : str, nonedited amino acid seq
    new_aa : str, edited amino acid seq
    start : start position of the first amino acid in the codon_window
    amino_acid_seq : dict, amino acid with {index : amino acid}

    Returns
    ------------
    mutation : str, translated formated mutation
    """
    result = []
    for i in range(len(aa)): 
        # if there is an edit made from base editing
        if aa[i] != new_aa[i]: 
            mut = aa[i] + str(start+i) + new_aa[i]
            # checks mutation against the protein sequence
            assert amino_acid_seq[start+i] == aa[i]
            # add edit to a list
            if mut not in result: 
                result.append(mut)
    mutation = '/'.join(result)
    return mutation

def categorize_mutations(mut_list): 
    """
    Categorizes mutations by Missense, Nonsense, Silent, etc
    
    Parameters
    ------------
    mut_list : list, a list of mutations the guide can access

    Returns
    ------------
    types : a list of unique mutations, sorted
    """
    types = []
    for mut in mut_list: 
        if len(mut) == 0 and 'Silent' not in types: 
            types.append('Silent')
        elif '.' in mut and 'Nonsense' not in types: 
            types.append('Nonsense')
        elif 'Missense' not in types: 
            types.append('Missense')
    return sorted(types)

def calc_target(row, window, mode, col_names): 
    """
    Finds the string of the bps in codons that may be edited
    Finds the string of the aa that may be edited

    Parameters
    ------------
    row : row of input df
    window : tuple, inclusive range for editing
    mode : str, 'DNA' or 'AA'
    col_names : tuple, (starting_frame, sgRNA_strand, gene_pos, seq_col)

    Returns
    ------------
    DNA or AA : str, DNA sequence of AA sequence that may be edited
    """
    frame = row[col_names[0]]
    guide = row[col_names[3]]
    num_aa = int(2+((window[1]-window[0]-1)//3))
    if row[col_names[1]] == 'sense': 
        start = (-1*frame)+3
        dna = guide[start:start+(num_aa*3)]
    else: 
        start = frame+1
        dna = rev_complement(complements, guide[start:start+(num_aa*3)])

    if mode == 'DNA': 
        return dna
    else: 
        return DNA_to_AA(dna, upper=False)
    
def calc_coding_window(row, window, col_names): 
    """
    Finds the string of the bps that may be edited
    
    Parameters
    ------------
    row : row of input df
    window : tuple, inclusive range for editing
    col_names : tuple, (starting_frame, sgRNA_strand, gene_pos, seq_col)

    Returns
    ------------
    seq : str, bps in codons that may be edited
    """
    if row[col_names[1]] == 'sense': 
        seq = row[col_names[3]][window[0]-1:window[1]]
    else: 
        seq = rev_complement(complements, row[col_names[3]][window[0]-1:window[1]])
    return seq

def calc_editing_window(row, window, col_names): 
    """
    Calculates the gene position of the window where an edit may occur
    
    Parameters
    ------------
    row : row of input df
    window : tuple, inclusive range for editing
    col_names : tuple, (starting_frame, sgRNA_strand, gene_pos, seq_col)

    Returns
    ------------
    gene_pos_window : tuple, gene position of window
    """
    if row[col_names[1]] == 'sense': 
        gene_pos_window = (row[col_names[2]]+window[0]-1, row[col_names[2]]+window[1]-1)
    else: 
        gene_pos_window = (row[col_names[2]]-window[0]+1, row[col_names[2]]-window[1]+1)
    return gene_pos_window
    