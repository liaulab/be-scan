"""
Author: Calvin XiaoYang Hu
Date: 231102

{Description: helper functions for processing guides and library data}
"""

import warnings
import re

from itertools import product
# from be_scan.sgrna._genomic_ import complements, rev_complement, DNA_to_AA
from _genomic_ import complements, rev_complement, DNA_to_AA

# FUNCTIONS FOR generate_library #

def filter_guide(g, PAM_regex, edit, window, excl_introns, excl_nonediting): 
    """
    Evaluates if a guide has a PAM and target residue within its window. 

    Parameters
    ------------
    excl_introns : bool, whether or not the editible base needs to be in an exon
    excl_nonediting : bool, whether or not the editible base needs to be present in the window
    """
    window_seq = g[0][window[0]-1:window[1]]
    if not excl_nonediting: 
        edit_in_window = True
    elif excl_introns: 
        edit_in_window = (edit[0] in window_seq)
    else: 
        edit_in_window = (edit[0].upper() in window_seq.upper())
    return (True if PAM_regex.match(g[1]) else False) and edit_in_window

def filter_sequence(results, sequence): 
    """
    Delete guides with TTTT which is a stop sequence for cloning
    """
    results = [r for r in results if not (sequence in r[0].upper())]
    return results

# FUNCTIONS FOR annotate #
def calc_target(row, window, col_names): 
    """
    Finds the string of the bps in codons that may be edited
    or Finds the string of the aa that may be edited
    """
    frame_col, strand_col, gene_pos_col, seq_col, window_start_col, window_end_col = col_names

    if row[window_start_col] == -1 and row[window_end_col] == -1: # NO CODING DNA IN WINDOW
        return None, None
    
    frame = row[frame_col]
    guide = row[seq_col]
    num_aa = int(2+((window[1]-window[0]-1)//3)) # pretty sure this works
    if row[strand_col] == 'sense': 
        start = (-1*frame)+3
        dna = guide[start:start+(num_aa*3)]
    else: 
        start = frame+1
        dna = rev_complement(complements, guide[start:start+(num_aa*3)])
    return dna, DNA_to_AA(dna, upper=False)
    
def calc_coding_window(row, window, col_names): 
    """
    Finds the fwd window of the bps that may be edited, whether it is exon or intron
    """
    frame_col, strand_col, gene_pos_col, seq_col, window_start_col, window_end_col = col_names

    if row[window_start_col] == -1 and row[window_end_col] == -1: # NO CODING DNA IN WINDOW
        return None, None

    frame = row[frame_col]
    if row[strand_col] == 'sense': start = (-1*frame)+3
    else:                          start = frame+1

    # RETURN EVEN IF THERE IS PARTIAL CODING DNA IN WINDOW #
    if row[strand_col] == 'sense': 
        return row[seq_col][window[0]-1:window[1]], int(window[0]-start-1)
    return rev_complement(complements, row[seq_col][window[0]-1:window[1]]), int(start-window[0])


def annotate_muts(row, edit, amino_acid_seq, col_names, pre, window, exons): 
    """
    Come up with list of annotations (ie F877L;F877P) for each guide
    """
    frame_col, strand_col, gene_pos_col, seq_col, window_start_col, window_end_col = col_names

    if row[window_start_col] == -1 and row[window_end_col] == -1: # NO CODING DNA IN WINDOW
        return None
    
    # EXTRACT DATA FROM ROW #
    frame, dir, pos, seq = row[frame_col], row[strand_col], row[gene_pos_col], row[seq_col]
    window1_pos, window2_pos = row[window_start_col], row[window_end_col]
    dna_window, dna, aa = row[f'{pre}_target_CDS'], row[f'{pre}_codon_window'], row[f'{pre}_residue_window']
    substitute_pos = int(row[f'{pre}_target_windowpos'])

    # USE POS TO CALCULATE WHICH AMINO ACID #
    if dir == 'sense' and pos != -1:               start = int((pos-frame)/3)+2
    elif dir == 'antisense' and pos != -1:         start = int((pos-frame+1)//3)-2
    elif dir == 'sense' and window1_pos != -1:
        if window1_pos < 3:                        start = 1
        else:                                      start = int((window1_pos-(window[0]+1)-frame)/3)+3
    elif dir == 'antisense' and window1_pos != -1: start = int((window1_pos+(window[0]-1)-frame+1)/3)-2
    elif dir == 'sense' and window2_pos != -1: 
        if window2_pos == 0:                       start = -1
        elif window2_pos <= 3:                     start = 0
        else:                                      start = int((window2_pos-(window[1]+1)-frame)/3)+3
    elif dir == 'antisense' and window2_pos != -1: start = int((window2_pos+(window[1]-1)-frame+1)/3)-2

    # LIST OF MUTATIONS FOR EACH ROW #
    mutation_details = []
    for m in mutation_combos(dna_window, edit, dir): 
        # COMPARE MUTATED DNA/AA WITH OLD DNA/AA #
        dna_temp, aa_temp = dna, aa
        # REPLACE CAUSED AN ISSUE FOR REPEATED AMINO ACIDS AND REGIONS #
        if dir == 'sense': new_dna = dna[:substitute_pos] + m + dna[substitute_pos+len(dna_window):]
        else: new_dna = dna[:len(dna)+substitute_pos+1-len(dna_window)] + m + dna[len(dna)+substitute_pos+1:]
        new_aa = DNA_to_AA(new_dna, upper=False)

        if window_on_boundary(new_dna): 
            dna_temp = fill_in_dna_exons(dna_temp, row, exons)
            aa_temp = DNA_to_AA(dna_temp, upper=False)
            new_dna = fill_in_dna_exons(new_dna, row, exons)
            new_aa = DNA_to_AA(new_dna, upper=False)

        # FORMAT MUTATIONS AS LETTER-NUMBER-LETTER #
        if dna_temp != new_dna: 
            mutations = format_mutation(aa_temp, new_aa, start, amino_acid_seq, dna_temp, new_dna, seq)
            mutation_details.append(mutations)
    return ';'.join(set(filter(None, mutation_details)))

def annotate_dual_muts(row, edit, amino_acid_seq, col_names, pre, window, exons): 
    """
    Come up with list of annotations (ie F877L;F877P) for each guide for a dual editor
    """
    frame_col, strand_col, gene_pos_col, seq_col, window_start_col, window_end_col = col_names

    if row[window_start_col] == -1 and row[window_end_col] == -1: # NO CODING DNA IN WINDOW
        return None
    
    # EXTRACT DATA FROM ROW #
    frame, dir, pos, seq = row[frame_col], row[strand_col], row[gene_pos_col], row[seq_col]
    window1_pos, window2_pos = row[window_start_col], row[window_end_col]
    dna_window, dna, aa = row[f'{pre}_target_CDS'], row[f'{pre}_codon_window'], row[f'{pre}_residue_window']
    substitute_pos = int(row[f'{pre}_target_windowpos'])

    # USE POS TO CALCULATE WHICH AMINO ACID #
    if dir == 'sense' and pos != -1:               start = int((pos-frame)/3)+2
    elif dir == 'antisense' and pos != -1:         start = int((pos-frame+1)//3)-2
    elif dir == 'sense' and window1_pos != -1:
        if window1_pos < 3:                        start = 1
        else:                                      start = int((window1_pos-(window[0]+1)-frame)/3)+3
    elif dir == 'antisense' and window1_pos != -1: start = int((window1_pos+(window[0]-1)-frame+1)/3)-2
    elif dir == 'sense' and window2_pos != -1: 
        if window2_pos == 0:                       start = -1
        elif window2_pos <= 3:                     start = 0
        else:                                      start = int((window2_pos-(window[1]+1)-frame)/3)+3
    elif dir == 'antisense' and window2_pos != -1: start = int((window2_pos+(window[1]-1)-frame+1)/3)-2

    # LIST OF MUTATIONS FOR EACH ROW #
    mutation_details = []
    edit1 = edit[0][0], edit[1][0]
    edit2 = edit[0][1], edit[1][1]
    for m1 in mutation_combos(dna_window, edit1, dir): 
        for m2 in mutation_combos(m1, edit2, dir): 
            # COMPARE MUTATED DNA/AA WITH OLD DNA/AA #
            dna_temp, aa_temp = dna, aa
            # REPLACE CAUSED AN ISSUE FOR REPEATED AMINO ACIDS AND REGIONS #
            if dir == 'sense': new_dna = dna[:substitute_pos] + m2 + dna[substitute_pos+len(dna_window):]
            else: new_dna = dna[:len(dna)+substitute_pos+1-len(dna_window)] + m2 + dna[len(dna)+substitute_pos+1:]
            new_aa = DNA_to_AA(new_dna, upper=False)

            if window_on_boundary(new_dna): 
                dna_temp = fill_in_dna_exons(dna_temp, row, exons)
                aa_temp = DNA_to_AA(dna_temp, upper=False)
                new_dna = fill_in_dna_exons(new_dna, row, exons)
                new_aa = DNA_to_AA(new_dna, upper=False)

            # FORMAT MUTATIONS AS LETTER-NUMBER-LETTER #
            if dna_temp != new_dna: 
                mutations = format_mutation(aa_temp, new_aa, start, amino_acid_seq, dna_temp, new_dna, seq)
                mutation_details.append(mutations)
    return ';'.join(set(filter(None, mutation_details)))

def fill_in_dna_exons(new_dna, row, exons): 
    new_dna_codons = [new_dna[i:i + 3] for i in range(0, len(new_dna), 3)]
    exon_i = row['exon']
    for i, codon in enumerate(new_dna_codons): 
        if window_on_boundary(codon): 
            num_lower = sum(1 for c in codon if c.islower())
            if codon[0].islower(): # LOOK AT PREVIOUS EXON #
                exon_seq = exons[exon_i - 1]
                replacement_codon = exon_seq[-num_lower:] + codon[-(3-num_lower):]
                new_dna_codons[i] = replacement_codon
            if codon[-1].islower(): # LOOK AT NEXT EXON #
                exon_seq = exons[exon_i + 1]
                replacement_codon = codon[:(3-num_lower)] + exon_seq[:num_lower]
                new_dna_codons[i] = replacement_codon
    return ''.join(new_dna_codons)

def window_on_boundary(dna): 
    """
    Function to determine whether a DNA sequence contains upper and lower case.
    """
    has_upper, has_lower = False, False
    for char in dna:
        if char.isupper():   has_upper = True
        elif char.islower(): has_lower = True
        if has_upper and has_lower:
            return True
    return False


def mutation_combos(guide_window, edit, dir):
    """
    Generates a list of all possible mutations of the guide_window
    """
    if dir == 'antisense': edit = complements[edit[0]], complements[edit[1]]
    mutated = []
    # CONVERT INPUT STRING TO A LIST #
    seq = list(guide_window)
    # FIND INDICES OF MUTATABLE BP IN SEQ #
    indices = [i for i, c in enumerate(seq) if c in edit[0]]

    # GENERATE KEY LETTER COMBOS AND INSERT INTO LIST #
    for t in product(edit[0]+edit[1], repeat=len(indices)):
        for i, c in zip(indices, t):
            seq[i] = c
        mutated.append(''.join(seq))
    return mutated

def format_mutation(aa, new_aa, start, amino_acid_seq, dna, new_dna, guide): 
    """
    Translates and formats the mutation with unedited-position-edited
    """
    result = []
    for i in range(len(aa)): 
        if dna[i*3:(i*3)+3] == new_dna[i*3:(i*3)+3]: # NO DNA MUT
            continue
        if aa[i] == '_': # INSIDE INTRON
            continue
        mut = aa[i] + str(start+i) + new_aa[i]
        assert start+i > 0, f'Error {start}'
        # CHECK MUT AGAINST PROTEIN SEQ #
        if amino_acid_seq is not None: 
            ### assert amino_acid_seq[start+i] == aa[i], f"Error: guide {guide}"
            if amino_acid_seq[start+i] != aa[i]: 
                warnings.warn(f"Error: guide {guide}")
        # ADD EDIT TO LIST #
        if mut not in result: 
            result.append(mut)
    mutation = '/'.join(result)
    return mutation

def determine_mutations(row, col_names, pre): 
    """
    Determine mutations based on predicted edits ie only Missense, Nonsense, Silent
    """
    frame_col, strand_col, gene_pos_col, seq_col, window_start_col, window_end_col = col_names

    if row[window_start_col] == -1 and row[window_end_col] == -1: # NO CODING DNA IN WINDOW
        return None
    
    types = []
    for muts in row[f'{pre}_mutations'].split(';'): 
        type = []
        muts_list = muts.split('/')
        for mut in muts_list: 
            if len(mut) == 0: 
                continue
            elif mut[0] == mut[-1]: type.append('Silent')
            elif '.' == mut[-1]: type.append('Nonsense')
            else: type.append('Missense')
        types.append('/'.join(type))
    return ';'.join(types)

def categorize_mutations(row, pre, col_names, window, edit_from): 
    """
    Categorizes mutations by predicted mutations and context metadata. 
    Based on a priority of: 
    1. Nonsense, 2. Splice Site, 3. Missense, 4. Intron, 5. Silent, 6. UTR, 7. Flank, 8. No Mutation
    """
    frame_col, strand_col, gene_pos_col, seq_col, window_start_col, window_end_col = col_names
    frame, dir, pos, seq = row[frame_col], row[strand_col], row[gene_pos_col], row[seq_col]

    # 'Splice-acceptor' means guide edits the 'AG' at 3' of intron
    # 'Splice-donor' means guide edits the 'GT' at 5' of intron
    if is_splice_acceptor_guide(row, pre, dir, window): 
        return 'Splice-acceptor'
    elif is_splice_donor_guide(row, pre, dir, window): 
        return 'Splice-donor'
    elif row[f'{pre}_muttypes'] is not None and 'Nonsense' in row[f'{pre}_muttypes']: 
        return 'Nonsense'
    elif row[f'{pre}_muttypes'] is not None and 'Missense' in row[f'{pre}_muttypes']: 
        return 'Missense'
    elif row[f'{pre}_muttypes'] is not None and 'Silent' in row[f'{pre}_muttypes']: 
        return 'Silent'
    else: 
        grna_seq = row[f'sgRNA_seq']
        seq_window = grna_seq[window[0]-1:window[1]]
        if any(char.islower() for char in seq_window) and edit_from.lower() in seq_window: 
            if row['UTR'] == True: 
                return 'UTR'
            else: 
                return 'Intron'
        else: 
            return 'No Mutation'

def is_splice_acceptor_guide(row, pre, dir, window): 
    pre_from = pre.split('to')[0]
    seq = row['sgRNA_seq']
    edit_wind = seq[window[0]-1:window[1]]
    if edit_wind is None: 
        return False
    if dir == 'sense' and bool(re.match(r'^[a-z]+ag[A-Z]+$', seq)):
        if ('A' in pre_from) and bool(re.match(r'^[a-z]*ag[A-Z]*$', edit_wind)) and seq[window[1]].isupper(): 
            return True
        elif ('A' in pre_from) and bool(re.match(r'^[a-z]*a$', edit_wind)) and seq[window[1]+1].isupper(): 
            return True
        if ('G' in pre_from) and bool(re.match(r'^[a-z]*ag[A-Z]*$', edit_wind)) and seq[window[1]].isupper(): 
            return True
        elif ('G' in pre_from) and bool(re.match(r'^g[A-Z]*$', edit_wind)) and seq[window[1]+1].isupper(): 
            return True
    elif dir == 'antisense' and bool(re.match(r'^[A-Z]+ct[a-z]+$', seq)):
        if ('T' in pre_from) and bool(re.match(r'^[A-Z]*ct[a-z]*$', edit_wind)) and seq[window[1]].islower(): 
            return True
        elif ('T' in pre_from) and bool(re.match(r'^t[a-z]*$', edit_wind)) and seq[window[1]+1].islower(): 
            return True
        if ('C' in pre_from) and bool(re.match(r'^[A-Z]*ct[a-z]*$', edit_wind)) and seq[window[1]].islower(): 
            return True
        elif ('C' in pre_from) and bool(re.match(r'^[A-Z]*c$', edit_wind)) and seq[window[1]+1].islower(): 
            return True
    return False

def is_splice_donor_guide(row, pre, dir, window): 
    pre_from = pre.split('to')[0]
    seq = row['sgRNA_seq']
    edit_wind = seq[window[0]-1:window[1]]
    if edit_wind is None: 
        return False
    if dir == 'sense' and bool(re.match(r'^[A-Z]*gt[a-z]*$', seq)):
        if ('T' in pre_from) and bool(re.match(r'^[A-Z]*gt[a-z]*$', edit_wind)) and seq[window[1]].islower(): 
            return True
        elif ('T' in pre_from) and bool(re.match(r'^t[a-z]*$', edit_wind)) and seq[window[1]+1].islower(): 
            return True
        if ('G' in pre_from) and bool(re.match(r'^[A-Z]*gt[a-z]*$', edit_wind)) and seq[window[1]].islower(): 
            return True
        elif ('G' in pre_from) and bool(re.match(r'^[A-Z]*g$', edit_wind)) and seq[window[1]+1].islower(): 
            return True
    elif dir == 'antisense' and bool(re.match(r'^[a-z]*ac[A-Z]*$', seq)):
        if ('A' in pre_from) and bool(re.match(r'^[a-z]*ac[A-Z]*$', edit_wind)) and seq[window[1]].isupper(): 
            return True
        elif ('A' in pre_from) and bool(re.match(r'^[a-z]*a$', edit_wind)) and seq[window[1]+1].isupper(): 
            return True
        if ('C' in pre_from) and bool(re.match(r'^[a-z]*ac[A-Z]*$', edit_wind)) and seq[window[1]].isupper(): 
            return True
        elif ('C' in pre_from) and bool(re.match(r'^c[A-Z]*$', edit_wind)) and seq[window[1]+1].isupper(): 
            return True
    return False

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
    if coding.isupper(): return "Exon"
    elif coding.islower(): return "Intron"
    else: return "Exon/Intron"

def annotate_intron_exon_anti(x, window): 
    coding = x[20-window[1]: 20-window[0]+1] # window=(4,8) so x[12:17]
    if coding.isupper(): return "Exon"
    elif coding.islower(): return "Intron"
    else: return "Exon/Intron"
