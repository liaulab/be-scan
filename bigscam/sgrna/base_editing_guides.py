"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: main functions for taking in a gene object and filtering potential guides based on 
              PAM, availability of target residue in editing window, and annotating the guides}
"""

import sys
import pandas as pd

sys.path.append('../helper')
from genomic import DNA_AA_map, base_editing_key, bases, complements, cas_key
from genomic import rev_complement, complement, protein_to_AAseq, process_PAM # DNA_to_AA

def find_gRNAs(aa_file, gene_object, mode, cas_type, target_codons=[], window=[4,8], PAM=None): 
    # mode: can be CBE or ABE, or just a list of 2 base edits ex ["C", "T"]
    # cas_type: Sp, SpG, SpRY
    # target_codons: list of codons that we want to make with our base edit
    # window: 4th to 8th bases inclusive by default, can be changed
    # PAM: optional field to input a custom PAM
    # Returns: a df of exon #, guides (23 bps), target (20 bps), fwd or rev, codon #, edit made
    amino_acid_seq = protein_to_AAseq(aa_file)
    
    # process mode
    if len(mode) == 3: 
        mode = base_editing_key[mode]
    assert len(mode) == 2
    
    # process PAM
    if PAM is None: 
        PAM = cas_key[cas_type]
    PAM_regex = process_PAM(PAM)
    
    # filter for PAM and contains editable base in window
    #    (seq, frame012 of first base, index of first base, exon number)
    fwd_results = [g.copy() for g in gene_object.fwd_guides if PAM_regex.match(g[0][-len(PAM):]) and 
                                                    mode[0] in g[0][window[0]-1:window[1]]]
    for g in fwd_results: 
        # mutates all residues according to the mode
        original = g[0][:12]
        mutated = g[0][:window[0]-1] + g[0][window[0]-1:window[1]].replace(mode[0], 
                                                                            mode[1]) + g[0][window[1]:12]
        assert(len(original)==len(mutated))
        # compares the residues to find which amino acids were altered and catalogs them
        start = (-1*g[1])+3
        original, mutated = original[start:start+9], mutated[start:start+9]
        original_aa, mutated_aa = '', ''
        edit, edit_ind = [], []
        for i in range(3): 
            if not original[i*3:(i+1)*3].isupper(): 
                continue
            original_aa += DNA_AA_map[original[i*3:(i+1)*3]]
            mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]

            assert amino_acid_seq[int((g[2]+1+start+(i*3))/3)+1]==original_aa[-1] # check we referenced correct aa
            if original_aa[-1] != mutated_aa[-1]: 
                edit.append(original_aa[-1] + ">" + mutated_aa[-1])
                edit_ind.append(int((g[2]+1+start+(i*3))/3)+1) # +1 for indexing starting at 1
                
        if len(edit) == 0: 
            edit.append('No Change')
        # append all information to dataframe
        g.append(original_aa)
        g.append(mutated_aa)
        g.append(edit)
        g.append(edit_ind)
        g.append('fwd')
        assert(len(g)) == 9
        
    # filter for PAM and contains editable base in window 
    #    (seq, frame012 of last base, index of last base, exon number)
    rev_results = [g.copy() for g in gene_object.rev_guides if PAM_regex.match(g[0][-len(PAM):]) and 
                                                    mode[0] in g[0][window[0]-1:window[1]]]
    for g in rev_results: 
        # mutates all residues according to the mode
        original = g[0][:12]
        mutated = g[0][:window[0]-1] + g[0][window[0]-1:window[1]].replace(mode[0], 
                                                                            mode[1]) + g[0][window[1]:12]
        assert(len(original)==len(mutated))
        # compares the residues to find which amino acids were altered and catalogs them
        original = rev_complement(complements, original[g[1]+1:g[1]+10])
        mutated = rev_complement(complements, mutated[g[1]+1:g[1]+10])
        original_aa, mutated_aa = '', ''
        edit, edit_ind = [], []
        for i in range(3): 
            if not original[i*3:(i+1)*3].isupper(): 
                continue
            original_aa += DNA_AA_map[original[i*3:(i+1)*3]]
            mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]
            assert amino_acid_seq[int((g[2]-10+1+(i*3))/3)+1]==original_aa[-1] # check we referenced correct aa
            if original_aa[-1] != mutated_aa[-1]: 
                edit.append(original_aa[-1] + ">" + mutated_aa[-1])
                edit_ind.append(int((g[2]-10+1+(i*3))/3))
                
        if len(edit) == 0: 
            edit.append('No Change')
        # append all information to dataframe
        g.append(original_aa)
        g.append(mutated_aa)
        g.append(edit)
        g.append(edit_ind)
        g.append('rev')
        assert(len(g)) == 9
        
    results = fwd_results + rev_results
    return pd.DataFrame(results)

