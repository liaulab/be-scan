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
    return PAM_regex.match(g[0][-len(PAM):]) and edit[0] in g[0][window[0]-1:window[1]]

# delete duplicate guides since duplicate guides are difficult to deconvolute
# input: a list of guides
# output a shorter list of guides with no duplicates
def filter_repeats(results): 
    seqs = [g[0] for g in results]
    dup_seqs = set([x for x in seqs if seqs.count(x) > 1])
    results = [g.copy() for g in results if g[0] not in dup_seqs]
    return results

def annotate_mutations(row, edit_from, edit_to, amino_acid_seq): 
    # extract relevant data from dataframe
    dir = row['sgRNA_strand']
    dna_window = row['target_CDS']
    dna, aa = row['codon_window'], row['residue_window']
    frame, pos = row['starting_frame'], row['gene_pos']

    # starting index of amino acid is different for sense vs anti
    if dir == 'sense': 
        start = int((pos+(-1*frame))/3)+2 # this is indexed at 1, for index at 0 +3 instead of +6
    else: 
        start = int((pos+(-1*frame)+1)/3)-2

    # add to a list of mutations for each row
    mutation_details = []
    for m in mutation_combos(dna_window, edit_from, edit_to, dir): 
        if m == dna_window: 
            continue
        # cnstruct new dna and new aa sequence from one mutation
        new_dna = dna.replace(dna_window, m)
        new_aa = DNA_to_AA(new_dna, upper=False)
        # write out which mutations are changed
        mutations = format_mutation(aa, new_aa, start, amino_acid_seq)
        mutation_details.append(mutations)
    return mutation_details

def mutation_combos(guide_window, edit_from, edit_to, dir):
    if dir == 'antisense': 
        edit_from = complements[edit_from]
        edit_to = complements[edit_to]
    mutated = []
    # Convert input string into a list so we can easily substitute letters
    seq = list(guide_window)
    # Find indices of key letters in seq
    indices = [ i for i, c in enumerate(seq) if c in edit_from]

    # Generate key letter combinations & place them into the list
    for t in product(edit_from+edit_to, repeat=len(indices)):
        for i, c in zip(indices, t):
            seq[i] = c
        mutated.append(''.join(seq))
    return mutated

def format_mutation(aa, new_aa, start, amino_acid_seq): 
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
    return '/'.join(result)

def categorize_mutations(mut_list): 
    types = []
    for mut in mut_list: 
        if len(mut) == 0 and 'Silent' not in types: 
            types.append('Silent')
        elif '.' in mut and 'Nonsense' not in types: 
            types.append('Nonsense')
        elif 'Missense' not in types: 
            types.append('Missense')
    return sorted(types)

def calc_target(row, window, mode): 
    frame = row['starting_frame']
    guide = row['sgRNA_seq']
    num_aa = int(2+((window[1]-window[0]-1)//3))
    if row['sgRNA_strand'] == 'sense': 
        start = (-1*frame)+3
        dna = guide[start:start+(num_aa*3)]
    else: 
        start = frame+1
        dna = rev_complement(complements, guide[start:start+(num_aa*3)])

    if mode == 'DNA': 
        return dna
    else: 
        return DNA_to_AA(dna, upper=False)
    
def calc_coding_window(x, window): 
    if x['sgRNA_strand'] == 'sense': 
        return x['sgRNA_seq'][window[0]-1:window[1]]
    else: 
        return rev_complement(complements, x['sgRNA_seq'][window[0]-1:window[1]])

def calc_editing_window(row, window): 
    if row['sgRNA_strand'] == 'sense': 
        return (row['gene_pos']+window[0]-1, row['gene_pos']+window[1]-1)
    else: 
        return (row['gene_pos']-window[0]+1, row['gene_pos']-window[1]+1)
    