"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: main functions for taking in a gene object and filtering potential guides based on 
              PAM, availability of target residue in editing window, and annotating the guides}
"""

import pandas as pd

from ._genomic_ import bases, complements, cas_key
from ._genomic_ import rev_complement, complement, protein_to_AAseq, process_PAM, make_mutations
from ._guides_ import filter_guide, filter_repeats
from ._aminoacid_ import find_aa_edits_fwd, find_aa_edits_rev


# this is the main function for taking in a gene object with all possible guides,
# then filtering based on given criteria
def identify_BE_guides(gene_object, cas_type, edit_from, edit_to, PAM=None, window=[4,8]): 
    # Parameters
    #    cas_type: Sp, SpG, SpRY, etc
    #    edit_from: the base (ACTG) to be replaced
    #    edit_to: the base (ACTG) to replace with
    #    window: editing window, 4th to 8th bases inclusive by default
    #    PAM: optional field to input a custom PAM
    #    Returns: a df of exon #, guides (23 bps), target (20 bps), fwd or rev
    # outputs: a list of fwd_guides and a list of rev_guides and the edit ie ['C', 'G']

    # process cas_type
    if cas_type not in list(cas_key.keys()): 
        raise Exception('Improper cas type input, the options are '+str(list(cas_key.keys())))
    
    assert edit_from in bases and edit_to in bases
    edit = edit_from, edit_to
    
    # process PAM, PAM input overrides cas_type
    if PAM is None: 
        PAM = cas_key[cas_type]
    PAM_regex = process_PAM(PAM)

    # process window
    assert window[1] >= window[0] and window[0] >= 0 
    assert window[1] <= len(gene_object.fwd_guides[0][0])

    # filter for PAM and contains editable base in window
    #    (seq, frame012 of first base, index of first base, exon number)
    fwd_results = [g.copy() for g in gene_object.fwd_guides if filter_guide(g, PAM_regex, PAM, edit, window)]

    # filter for PAM and contains editable base in window 
    #    (seq, frame012 of last base, index of last base, exon number)
    rev_results = [g.copy() for g in gene_object.rev_guides if filter_guide(g, PAM_regex, PAM, edit, window)]

    # filter out repeating guides in fwd_results list
    fwd_results = filter_repeats(fwd_results)
    # filter out repeating guides in rev_results list
    rev_results = filter_repeats(rev_results)

    ### filter out guides with off target editing in the genome


    return fwd_results, rev_results, edit

# this is the main function for taking in lists of guides, 
# then annotating all their predicted edits
def annotate_BE_guides(protein_filepath, fwd_guides, rev_guides, edit_from, edit_to, window=[4,8]): 
    ### ADD DOCS
    # Parameters
    #    protein_filepath: filepath to an amino acid sequence corresponding to gene file
    #    fwd_guides, rev_guides: generated from identify_guides
    #    edit_from: the base (ACTG) to be replaced
    #    edit_to: the base (ACTG) to replace with
    #    window: editing window, 4th to 8th bases inclusive by default
    # outputs: a dataframe
    
    ### target_codons: list of codons that we want to make with our base edit

    # codon indices, predicted edits made
    amino_acid_seq = protein_to_AAseq(protein_filepath)
    num_aa = 3 # this is the num of amino acids we look ahead in our frame

    for g in fwd_guides: 
        # mutates all residues according to the mode, every combination of residue mutations
        original = g[0][:12] # a string
        guide_window = g[0][window[0]-1:window[1]] # a string
        mutateds = [g[0][:window[0]-1] + m + g[0][window[1]:12] for m in make_mutations(guide_window, edit_from, edit_to)] # list of strings 

        # compares the residues to find which amino acids were altered and catalogs them
        edits, edit_inds = [], [] # lists of lists, of all edits for all possible mutations
        start = (-1*g[1])+3
        orig = original[start:start+(num_aa*3)]
        for m in mutateds: 
            # for each possible mutation, come up with the list of amino acid changes
            edit, edit_ind = find_aa_edits_fwd(m, g, start, orig, num_aa, amino_acid_seq)
            edits.append(edit)
            edit_inds.append(edit_ind)
        # append all information to dataframe
        g.extend([edits, edit_inds, 'fwd'])
        assert(len(g)) == 7
        
    for g in rev_guides: 
        # mutates all residues according to the mode, every combination of residue mutations
        original = g[0][:12] # a string
        guide_window = g[0][window[0]-1:window[1]] # a string
        mutateds = [g[0][:window[0]-1] + m + g[0][window[1]:12] for m in make_mutations(guide_window, edit_from, edit_to)] # a list of strings

        # compares the residues to find which amino acids were altered and catalogs them
        edits, edit_inds = [], [] # lists of lists, of all edits for all possible mutations
        start = g[1]+1
        orig = rev_complement(complements, original[start:start+(num_aa*3)])
        for m in mutateds: 
            # for each possible mutation, come up with the list of amino acid changes
            edit, edit_ind = find_aa_edits_rev(m, g, start, orig, num_aa, amino_acid_seq)
            edits.append(edit)
            edit_inds.append(edit_ind)
        # append all information to dataframe
        g.extend([edits, edit_inds, 'rev'])
        assert(len(g)) == 7
        
#     print(pd.DataFrame(fwd_guides + rev_guides))
    return pd.DataFrame(fwd_guides + rev_guides)
