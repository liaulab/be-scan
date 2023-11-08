"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: main functions for taking in a gene object and filtering potential guides based on 
              PAM, availability of target residue in editing window, and annotating the guides}
"""

import pandas as pd

from ._genomic_ import DNA_AA_map, base_editing_key, bases, complements, cas_key
from ._genomic_ import rev_complement, complement, protein_to_AAseq, process_PAM # DNA_to_AA

def identify_guides(gene_object, cas_type, mode, PAM=None, window=[4,8]): 
    # Parameters
    #    mode: can be CBE or ABE
    #    cas_type: Sp, SpG, SpRY, etc
    #    window: 4th to 8th bases inclusive by default, can be changed
    #    PAM: optional field to input a custom PAM
    #    Returns: a df of exon #, guides (23 bps), target (20 bps), fwd or rev

    # process cas_type
    if cas_type not in list(cas_key.keys()): 
        raise Exception('Improper cas type input, the options are '+str(list(cas_key.keys())))
    
    # process mode
    if mode not in list(base_editing_key.keys()): 
        raise Exception('Improper mode input, the options are '+str(list(base_editing_key.keys())))
    edit = base_editing_key[mode]    
    
    # process PAM, PAM input overrides cas_type
    if PAM is None: 
        PAM = cas_key[cas_type]
    PAM_regex = process_PAM(PAM)

    # process window
    assert window[1] >= window[0] and window[0] >= 0 and window[1] <= len(gene_object.fwd_guides[0][0])

    # filter for PAM and contains editable base in window
    #    (seq, frame012 of first base, index of first base, exon number)
    fwd_results = [g.copy() for g in gene_object.fwd_guides if PAM_regex.match(g[0][-len(PAM):]) and 
                                                    edit[0] in g[0][window[0]-1:window[1]]]

    # filter for PAM and contains editable base in window 
    #    (seq, frame012 of last base, index of last base, exon number)
    rev_results = [g.copy() for g in gene_object.rev_guides if PAM_regex.match(g[0][-len(PAM):]) and 
                                                    edit[0] in g[0][window[0]-1:window[1]]]

    return fwd_results, rev_results, edit


def annotate_guides(protein_filepath, fwd_guides, rev_guides, mode, window=[4,8]): 
    # target_codons: list of codons that we want to make with our base edit

    # codon indices, predicted edits made
    amino_acid_seq = protein_to_AAseq(protein_filepath)

    for g in fwd_guides: 
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
        
    for g in rev_guides: 
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
        
    results = fwd_guides + rev_guides
    return pd.DataFrame(results)

