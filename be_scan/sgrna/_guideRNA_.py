"""
Author: Calvin XiaoYang Hu
Date: 231102

{Description: helper functions for processing guides and library data}
"""

import re
from itertools import product
from be_scan.sgrna._genomic_ import complements, rev_complement, DNA_to_AA
from be_scan.sgrna._gene_ import GeneForCRISPR

# evaluates if guide has PAM and has a residue in window
# returns TRUE or FALSE
def filter_guide(g, PAM_regex, edit, window, exclude_introns, exclude_nontargeting): 
    """
    Evaluates if a guide has a PAM and target residue within its window. 

    Parameters
    ------------
    g : str, guide sequence
    PAM_regex : regex str, regez string constructued from PAM
    PAM : str, PAM sequence
    edit : tuple, (edit_from edit_to)
    window : tuple, inclusive range for editing
    exclude_introns : bool, whether or not the editible base needs to be in an intron
    exclude_nontargeting : bool, whether or not the editible base needs to be in the window

    Returns
    ------------
    bool True or False
    """
    seq = g[0]
    window = seq[window[0]-1:window[1]]
    if not exclude_nontargeting: 
        edit_in_window = True
    elif exclude_introns: 
        edit_in_window = (edit[0] in window)
    else: 
        edit_in_window = (edit[0].upper() in window) or (edit[0].lower() in window)
    return (True if PAM_regex.match(g[1]) else False) and edit_in_window

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

def annotate_mutations(row, edit, amino_acid_seq, col_names, prefix): 
    """
    Come up with list of annotations (ie F877L, F877P, F877L/F877P)
    for each guide. 

    Parameters
    ------------
    row : row of input df
    edit : tuple, (edit_from edit_to) each are one character
    amino_acid_seq : dict, amino acid with {index : amino acid}
    col_names : tuple, (starting_frame, sgRNA_strand, gene_pos, seq_col)

    Returns
    ------------
    mutation_details : list, a list of mutations
    """
    assert len(edit[0]) == 1 and len(edit[1]) == 1
    # if guide position is unannotated
    if row[col_names[2]] == -1: 
        return None
    
    # extract relevant data from dataframe
    frame, dir, pos = row[col_names[0]], row[col_names[1]], row[col_names[2]]
    dna_window, dna, aa = row[prefix+'_target_CDS'], row[prefix+'_codon_window'], row[prefix+'_residue_window']

    # starting index of amino acid is different for sense vs anti
    if dir == 'sense': 
        start = int((pos+(-1*frame))/3)+2 # indexed at 1, for index at 0 +3 instead of +6
    else: 
        start = int((pos+(-1*frame)+1)/3)-2

    # add to a list of mutations for each row
    mutation_details = []
    for m in mutation_combos(dna_window, edit, dir): 
        # construct new dna and new aa sequence from one potential mutation
        new_dna = dna.replace(dna_window, m)
        new_aa = DNA_to_AA(new_dna, upper=False)
        # write out which mutations are changed
        mutations = format_mutation(aa, new_aa, start, amino_acid_seq, row[col_names[3]])
        mutation_details.append(mutations)
    return mutation_details[1:]

### only works for dual editor for now
def annotate_dual_mutations(row, edit, amino_acid_seq, col_names, prefix): 
    """
    Come up with list of annotations (ie F877L, F877P, F877L/F877P)
    for each guide for multiple possible bp edits. 

    Parameters
    ------------
    row : row of input df
    edit : tuple, (edit_from edit_to) each are multiple characters
    amino_acid_seq : dict, amino acid with {index : amino acid}
    col_names : tuple, (starting_frame, sgRNA_strand, gene_pos, seq_col)

    Returns
    ------------
    mutation_details : list, a list of mutations
    """
    assert len(edit[0]) > 1 and len(edit[1]) > 1
    # if guide position is unannotated
    if row[col_names[2]] == -1: 
        return None
    
    # extract relevant data from dataframe
    frame, dir, pos = row[col_names[0]], row[col_names[1]], row[col_names[2]]
    dna_window, dna, aa = row[prefix+'_target_CDS'], row[prefix+'_codon_window'], row[prefix+'_residue_window']

    # starting index of amino acid is different for sense vs anti
    if dir == 'sense': 
        start = int((pos+(-1*frame))/3)+2 # indexed at 1, for index at 0 +3 instead of +6
    else: 
        start = int((pos+(-1*frame)+1)/3)-2

    # add to a list of mutations for each row
    mutation_details = []
    edit1 = (edit[0][0], edit[1][0])
    edit2 = (edit[0][1], edit[1][1])
    for m1 in mutation_combos(dna_window, edit1, dir): 
        for m2 in mutation_combos(m1, edit2, dir): 
            # construct new dna and new aa sequence from one potential mutation
            new_dna = dna.replace(dna_window, m2)
            new_aa = DNA_to_AA(new_dna, upper=False)
            # write out which mutations are changed
            mutations = format_mutation(aa, new_aa, start, amino_acid_seq, row[col_names[3]])
            mutation_details.append(mutations)
    return mutation_details[1:]

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

def format_mutation(aa, new_aa, start, amino_acid_seq, x): 
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
            if amino_acid_seq is not None: 
                assert amino_acid_seq[start+i] == aa[i], 'guides '+f"{x}" + ", " + aa[i] + str(start+i) + new_aa[i] + ", " + amino_acid_seq[start+i]
            # add edit to a list
            if mut not in result: 
                result.append(mut)
    mutation = '/'.join(result)
    return mutation

def categorize_mutations(row, col_names, prefix): 
    """
    Categorizes mutations by Missense, Nonsense, Silent, etc
    
    Parameters
    ------------
    mut_list : list, a list of mutations the guide can access

    Returns
    ------------
    types : a list of unique mutations, sorted
    """
    if row[col_names[2]] == -1: 
        return None
    
    types = []
    for mut in row[prefix+'_mutations']: 
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
    if row[col_names[2]] == -1: 
        return None
    
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
    if row[col_names[2]] == -1: 
        return None

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
    if row[col_names[2]] == -1: 
        return None
    
    if row[col_names[1]] == 'sense': 
        gene_pos_window = (row[col_names[2]]+window[0]-1, row[col_names[2]]+window[1]-1)
    else: 
        gene_pos_window = (row[col_names[2]]-window[0]+1, row[col_names[2]]-window[1]+1)
    return gene_pos_window
    
def parse_position_frames(guides_df, seq_col, gene_filepath, 
                          frame_col, strand_col, gene_pos_col): 
    """
    Populates input df with the start position of the guide, the 
    codon frame of the starting guide, and the strand direction. 
    This information is necessary for downstream annotation. 
    
    Parameters
    ------------
    guides_df : pandas df, only needs to contain a column seq_col
    seq_col : str, name of column with sgRNA sequence
    gene_filepath : str or path, filepath to the gene
    frame_col : str, name of column with codon frame
    strand_col : str, name of column with direction
    gene_pos_col : str, name of column with position on the gene
    
    Returns
    ------------
    guides_df : populated with seq_col, frame_col, strand_col, gene_pos_col
    """
    # extract guides, this is the only required column
    guides = guides_df[seq_col]

    # make gene object to string match against
    gene = GeneForCRISPR(filepath=gene_filepath)
    gene.parse_exons()
    gene_seq = ''.join(gene.exons)

    frames, positions, strands = [], [], []
    for guide in guides: 
        # find if coding region fwd or rev matches any part of the gene
        coding_guide_fwd = extract_uppercase_letters(guide)
        coding_guide_rev = rev_complement(complements, coding_guide_fwd)
        pos_fwd = list(re.finditer(coding_guide_fwd, gene_seq))
        pos_rev = list(re.finditer(coding_guide_rev, gene_seq))
        fwd_ind = find_first_uppercase_index(guide) # only need fwd bc rev_complement was taken

        # guide only matches one position in sense OR antisense (ideal)
        if len(pos_fwd) == 1 and len(pos_rev) == 0: 
            pos = pos_fwd[0].start()-fwd_ind
            stran = 'sense'
            frame = pos%3
        elif len(pos_fwd) == 0 and len(pos_rev) == 1: 
            pos = pos_rev[0].end()+fwd_ind-1
            stran = 'antisense'
            frame = (pos)%3
        else: 
            # guide matches sequences sense AND antisense, use first sense match
            if len(pos_fwd) > 0 and len(pos_rev) > 0: 
                print('The guide', guide, 'match sense and antisense strands.')
            # guide matches multiple positions in sense OR antisense, use first match from either
            elif len(pos_fwd) > 1 or len(pos_rev) > 1: 
                print('The guide', guide, 'has many occurrences.')
            else: 
                print('The guide', guide, 'couldn\'t be matched.')
            pos = -1
            stran = 'unknown'
            frame = -1

        # add information to lists
        positions.append(pos)
        strands.append(stran)
        frames.append(frame)

    # set new columns and return df
    if gene_pos_col not in guides_df.columns: 
        guides_df[gene_pos_col] = positions
    if frame_col not in guides_df.columns: 
        guides_df[frame_col] = frames
    if strand_col not in guides_df.columns: 
        guides_df[strand_col] = strands
    return guides_df

def extract_uppercase_letters(input):
    """
    Returns uppercase letters of input string.
    """
    return ''.join(char for char in input if char.isupper())

def find_first_uppercase_index(input_string):
    """
    Finds index of the first uppercase letter. 
    """
    for index, char in enumerate(input_string):
        if char.isupper():
            return index
        
def annotate_intron_exon(x, window): 
    coding = x[window[0]-1: window[1]]
    if coding.isupper(): 
        return "Exon"
    elif coding.islower(): 
        return "Intron"
    else: 
        return "Exon/Intron"