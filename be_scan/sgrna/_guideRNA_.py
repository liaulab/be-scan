"""
Author: Calvin XiaoYang Hu
Date: 231102

{Description: helper functions for processing guides and library data}
"""

import warnings
import re

from itertools import product
from be_scan.sgrna._genomic_ import complements, rev_complement, DNA_to_AA
# from _genomic_ import complements, rev_complement, DNA_to_AA

# FUNCTIONS FOR generate_library #

# All guides apart from exonic ones: set to –1
# Exonic guides: assign position based on 1-indexed protein sequence
# No C: Use the center of the editing window (convert to residue number)
# Yes C:
    # Silent: Use the center of the editing window
    # Missense: Take mutated residue numbers and average
    # Nonsense: Take nonsense mutation residue numbers (avg if multiple)
def assign_position(mutations, start_pos, end_pos, gene, muttype):
    if gene == 'Control' or gene == 'Essential':
        return -1

    # By default, all -1
    if muttype in ['Splice-donor', 'Splice-acceptor', 'UTR', 'Intron']:
        return -1

    # for Nonsense, Missense, Silent
    if not isinstance(mutations, float) and mutations is not None:

        nonsense_list = []
        missense_list = []
        silent_list = []

        mutations_list = re.split("/|;", mutations)
        for mutation in mutations_list:
          if len(mutation) == 0:
              continue
          numerals = re.findall(r'\d+', mutation)
          num = int(''.join(numerals))
          if mutation[-1] == '.' and mutation[0] != '.':
              nonsense_list.append(num)
          elif mutation[0] == mutation[-1]:
              silent_list.append(num)
          else:
              missense_list.append(num)

        nonsense_list = list(set(nonsense_list))
        missense_list = list(set(missense_list))
        silent_list = list(set(silent_list))

        if len(nonsense_list) > 0:
            return float(sum(nonsense_list)/len(nonsense_list))
        elif len(missense_list) > 0:
            return float(sum(missense_list)/len(missense_list))
        elif len(silent_list) > 0:
            return float(sum(silent_list)/len(silent_list))

    # for Splice-acceptor, Splice-donor, UTR, Intron, No Mutation
    #    Splice only if Missense edits are not present
    #    Splice takes precedence over Missense etc for categorization, but not position labeling

    # For No Mutation, take center of editing window
    if start_pos != -1 and end_pos != -1:
        center = (start_pos + end_pos) / 2
        return 1+(center//3)
    # elif start_pos != -1:
    #   return 1+(start_pos)//3
    # elif end_pos != -1:
    #   return 1+(end_pos)//3
    else:
        return -1


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

def codon_window_geometry(frame, strand, window, guide_len):
    """
    Where the translated codon window sits inside a guide, and how long it is.

    The editing window is `window` (1-indexed, inclusive), so it occupies guide
    indices window[0]-1 .. window[1]-1. The codon window has to begin on a
    codon boundary at or before the editing window and extend past its end.

    Codon boundaries depend on the strand:
      sense     -- guide index k is a boundary when (k + frame) % 3 == 0
      antisense -- boundaries are offset by frame+1

    Both rules reproduce what the previous hardcoded expressions produced at
    the default window [4, 8]; they are simply written in terms of `window`
    rather than baked in as constants. num_aa is likewise held at its original
    value and only grown when a wider window needs more codons, so the
    geometry at [4, 8] is unchanged.

    Returns (start, num_aa): the guide index the codon window begins at, and
    how many codons it spans.
    """
    first = window[0] - 1          # guide index of the window's first base
    last_exclusive = window[1]     # one past the window's last base

    if strand == 'sense': phase = (-frame) % 3
    else:                 phase = (frame + 1) % 3

    # largest codon boundary at or before the start of the editing window
    start = first - ((first - phase) % 3)

    # a codon beginning before the guide cannot be translated from the guide
    # alone, so step forward to the first boundary inside it
    while start < 0:
        start += 3

    # keep the original width, widening only when the window demands it
    num_aa = int(2+((window[1]-window[0]-1)//3))
    needed = -(-(last_exclusive - start) // 3)   # ceiling division
    num_aa = max(num_aa, needed)

    # never read past the end of the guide
    while num_aa > 0 and start + (num_aa * 3) > guide_len:
        num_aa -= 1

    return start, num_aa


def first_residue_number(strand, frame, pos, window1_pos, window2_pos,
                         window, guide_len, dna_len):
    """
    Residue number of the first codon in the translated codon window.

    Everything downstream indexes off this: format_mutation emits
    aa[i] + str(start+i), so if this is off by one then every mutation call
    for the guide is misnumbered.

    A guide carries up to three coordinates, any of which is -1 when the
    corresponding base sits in an intron:
      pos          -- CDS index of guide base 0
      window1_pos  -- CDS index of the editing window's first base
      window2_pos  -- CDS index of the editing window's last base

    Each is converted to the CDS index of the codon window's first base, then
    to a 1-indexed residue number. On the antisense strand the CDS coordinate
    runs backwards relative to the guide, so the codon window's first base in
    reading order is its last base in guide order.

    Reproduces the previous hardcoded arithmetic exactly at window [4, 8].
    """
    start, _ = codon_window_geometry(frame, strand, window, guide_len)

    if strand == 'sense':
        # CDS index of guide base 0, from whichever coordinate is available
        if pos != -1:           base = pos
        elif window1_pos != -1: base = window1_pos - (window[0]-1)
        elif window2_pos != -1: base = window2_pos - (window[1]-1)
        else:                   return -1
        cds_index = base + start
    else:
        # antisense: guide base k sits at CDS index (base - k)
        if pos != -1:           base = pos
        elif window1_pos != -1: base = window1_pos + (window[0]-1)
        elif window2_pos != -1: base = window2_pos + (window[1]-1)
        else:                   return -1
        # first base in reading order is the far end of the codon window
        cds_index = base - (start + dna_len - 1)

    return (cds_index // 3) + 1


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
    start, num_aa = codon_window_geometry(frame, row[strand_col], window, len(guide))
    if num_aa == 0:
        return None, None
    if row[strand_col] == 'sense':
        dna = guide[start:start+(num_aa*3)]
    else:
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
    guide = row[seq_col]
    start, num_aa = codon_window_geometry(frame, row[strand_col], window, len(guide))
    if num_aa == 0:
        return None, None

    # CLIP THE EDITING WINDOW TO THE CODON WINDOW #
    # A base outside the translated codons cannot be substituted into them.
    # This only bites when the editing window reaches past either end of the
    # codon window; at the default [4, 8] the clip is a no-op and the values
    # below are identical to the unclipped ones.
    lo = max(window[0]-1, start)
    hi = min(window[1]-1, start+(num_aa*3)-1)
    if hi < lo:
        return None, None

    # RETURN EVEN IF THERE IS PARTIAL CODING DNA IN WINDOW #
    if row[strand_col] == 'sense':
        return guide[lo:hi+1], int(lo-start)
    return rev_complement(complements, guide[lo:hi+1]), int(start-lo-1)


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
    start = first_residue_number(
        dir, frame, pos, window1_pos, window2_pos, window, len(seq), len(dna))

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
    start = first_residue_number(
        dir, frame, pos, window1_pos, window2_pos, window, len(seq), len(dna))

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
            # BOUNDS-CHECKED: exons[exon_i-1] would wrap to the last exon and
            # silently rebuild the codon from the wrong end of the gene, and
            # exons[exon_i+1] would raise past the final exon
            if codon[0].islower() and exon_i-1 >= 0: # LOOK AT PREVIOUS EXON #
                exon_seq = exons[exon_i - 1]
                replacement_codon = exon_seq[-num_lower:] + codon[-(3-num_lower):]
                new_dna_codons[i] = replacement_codon
            if codon[-1].islower() and exon_i+1 < len(exons): # LOOK AT NEXT EXON #
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

def _is_upper_at(seq, index):
    """
    Whether seq[index] is uppercase, False when the index is off either end.

    The splice detectors look one and two bases past the editing window to
    confirm which side of the junction the window sits on. Those look-aheads
    run off the end of the guide once the window reaches its 3' end, so the
    access is bounds-checked. Out of range means there is no evidence of a
    junction, which is the same answer as "not a splice guide".
    """
    return 0 <= index < len(seq) and seq[index].isupper()


def _is_lower_at(seq, index):
    """Whether seq[index] is lowercase, False when the index is off either end."""
    return 0 <= index < len(seq) and seq[index].islower()


def is_splice_acceptor_guide(row, pre, dir, window):
    pre_from = pre.split('to')[0]
    seq = row['sgRNA_seq']
    edit_wind = seq[window[0]-1:window[1]]
    if edit_wind is None: 
        return False
    if dir == 'sense' and bool(re.match(r'^[a-z]+ag[A-Z]+$', seq)):
        if ('A' in pre_from) and bool(re.match(r'^[a-z]*ag[A-Z]*$', edit_wind)) and _is_upper_at(seq, window[1]): 
            return True
        elif ('A' in pre_from) and bool(re.match(r'^[a-z]*a$', edit_wind)) and _is_upper_at(seq, window[1]+1): 
            return True
        if ('G' in pre_from) and bool(re.match(r'^[a-z]*ag[A-Z]*$', edit_wind)) and _is_upper_at(seq, window[1]): 
            return True
        elif ('G' in pre_from) and bool(re.match(r'^g[A-Z]*$', edit_wind)) and _is_upper_at(seq, window[1]+1): 
            return True
    elif dir == 'antisense' and bool(re.match(r'^[A-Z]+ct[a-z]+$', seq)):
        if ('T' in pre_from) and bool(re.match(r'^[A-Z]*ct[a-z]*$', edit_wind)) and _is_lower_at(seq, window[1]): 
            return True
        elif ('T' in pre_from) and bool(re.match(r'^t[a-z]*$', edit_wind)) and _is_lower_at(seq, window[1]+1): 
            return True
        if ('C' in pre_from) and bool(re.match(r'^[A-Z]*ct[a-z]*$', edit_wind)) and _is_lower_at(seq, window[1]): 
            return True
        elif ('C' in pre_from) and bool(re.match(r'^[A-Z]*c$', edit_wind)) and _is_lower_at(seq, window[1]+1): 
            return True
    return False

def is_splice_donor_guide(row, pre, dir, window): 
    pre_from = pre.split('to')[0]
    seq = row['sgRNA_seq']
    edit_wind = seq[window[0]-1:window[1]]
    if edit_wind is None: 
        return False
    if dir == 'sense' and bool(re.match(r'^[A-Z]*gt[a-z]*$', seq)):
        if ('T' in pre_from) and bool(re.match(r'^[A-Z]*gt[a-z]*$', edit_wind)) and _is_lower_at(seq, window[1]): 
            return True
        elif ('T' in pre_from) and bool(re.match(r'^t[a-z]*$', edit_wind)) and _is_lower_at(seq, window[1]+1): 
            return True
        if ('G' in pre_from) and bool(re.match(r'^[A-Z]*gt[a-z]*$', edit_wind)) and _is_lower_at(seq, window[1]): 
            return True
        elif ('G' in pre_from) and bool(re.match(r'^[A-Z]*g$', edit_wind)) and _is_lower_at(seq, window[1]+1): 
            return True
    elif dir == 'antisense' and bool(re.match(r'^[a-z]*ac[A-Z]*$', seq)):
        if ('A' in pre_from) and bool(re.match(r'^[a-z]*ac[A-Z]*$', edit_wind)) and _is_upper_at(seq, window[1]): 
            return True
        elif ('A' in pre_from) and bool(re.match(r'^[a-z]*a$', edit_wind)) and _is_upper_at(seq, window[1]+1): 
            return True
        if ('C' in pre_from) and bool(re.match(r'^[a-z]*ac[A-Z]*$', edit_wind)) and _is_upper_at(seq, window[1]): 
            return True
        elif ('C' in pre_from) and bool(re.match(r'^c[A-Z]*$', edit_wind)) and _is_upper_at(seq, window[1]+1): 
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

def annotate_intron_exon_anti(x, window, guide_len=None):
    """
    Where the mirrored editing window sits on an antisense guide.

    The window is mirrored because an antisense guide reads in the opposite
    direction to the coding sequence it is compared against. guide_len
    defaults to len(x), which for a standard 20 nt spacer gives x[12:17] at
    window (4,8) -- identical to the length this previously hardcoded.
    """
    if guide_len is None:
        guide_len = len(x)
    coding = x[guide_len-window[1]: guide_len-window[0]+1] # window=(4,8) so x[12:17]
    if coding.isupper(): return "Exon"
    elif coding.islower(): return "Intron"
    else: return "Exon/Intron"
