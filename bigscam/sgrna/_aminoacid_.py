"""
Author: Calvin XiaoYang Hu
Date: 231103

{Description: some base pair to amino acid translation functions}
"""

from ._genomic_ import rev_complement, complements

# translating DNA sequence to amino acid sequence
DNA_AA_map = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
              "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
              "TAT":"Y", "TAC":"Y", "TAA":".", "TAG":".",
              "TGT":"C", "TGC":"C", "TGA":".", "TGG":"W",
              "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
              "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
              "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
              "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
              "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
              "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
              "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
              "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
              "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
              "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
              "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
              "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G", }

# function to find all aa edit in the fwd direction
# given 1 mutant and information on where to start translating
def find_aa_edits_fwd(m, g, start, orig, num_aa, amino_acid_seq): 
    mutated = m[start:start+(num_aa*3)] # the mutated string that we are converting
    original_aa, mutated_aa = '', ''
    edit, edit_ind = [], []
    # for the next num_aa amino acids, log each translation and log the changes if there are any
    for i in range(num_aa): 
        if not orig[i*3:(i+1)*3].isupper(): # make sure bps are all coding (caps)
            continue
        # add translations to a string
        original_aa += DNA_AA_map[orig[i*3:(i+1)*3]]
        mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]
        assert amino_acid_seq[int((g[2]+1+start+(i*3))/3)+1]==original_aa[-1] # check we referenced correct aa
        # log all changed aa
        if original_aa[-1] != mutated_aa[-1]: 
            edit.append(original_aa[-1] + ">" + mutated_aa[-1])
            edit_ind.append(int((g[2]+1+start+(i*3))/3)+1) # +1 for indexing starting at 1
            
    if len(edit) == 0: 
        edit.append('No Change')
    return edit, edit_ind
    
# function to find all aa edit in the fwd direction
# given 1 mutant and information on where to start translating
def find_aa_edits_rev(m, g, start, orig, num_aa, amino_acid_seq): 
    mutated = rev_complement(complements, m[start:start+(num_aa*3)])
    original_aa, mutated_aa = '', ''
    edit, edit_ind = [], []
    # for the next num_aa amino acids, log each translation and log the changes if there are any
    for i in range(num_aa): 
        if not orig[i*3:(i+1)*3].isupper(): # make sure bps are all coding (caps)
            continue
        # add translations to a string
        original_aa += DNA_AA_map[orig[i*3:(i+1)*3]]
        mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]
        assert amino_acid_seq[int((g[2]-10+1+(i*3))/3)+1]==original_aa[-1] # check we referenced correct aa
        # log all changed aa
        if original_aa[-1] != mutated_aa[-1]: 
            edit.append(original_aa[-1] + ">" + mutated_aa[-1])
            edit_ind.append(int((g[2]-10+1+(i*3))/3))
            
    if len(edit) == 0: 
        edit.append('No Change')
    return edit, edit_ind
