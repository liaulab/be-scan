#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 15 12:52:34 2021

@author: nlclarinet8, ceejaylee

Description: This .py file contains the main commands to run to generate CBE analysis. This was written mostly by Nick Lue, with modifications by Ceejay Lee to integrate better with the lab's gRNA-generating pipeline.
"""

# Import Packages

import os
import pathlib
import time
import sys
import csv
from collections import OrderedDict, Counter
import itertools
import pdb
import math

import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from ._utils import (
    processIFile, find_Indexes, generateSubsets, generateMutCom, get_aa, find_sgRNA_pos, rev_complement, convert_guideInput
)

# this is the main function for taking in lists of guides, 
# then annotating all their predicted edits
# def annotate_BE_guides(protein_filepath, fwd_guides, rev_guides, edit_from, edit_to, window=[4,8]): 
#     ### ADD DOCS
#     # Parameters
#     #    protein_filepath: filepath to an amino acid sequence corresponding to gene file
#     #    fwd_guides, rev_guides: generated from identify_guides
#     #    edit_from: the base (ACTG) to be replaced
#     #    edit_to: the base (ACTG) to replace with
#     #    window: editing window, 4th to 8th bases inclusive by default
#     # outputs: a dataframe
    
#     ### target_codons: list of codons that we want to make with our base edit

#     # codon indices, predicted edits made
#     amino_acid_seq = protein_to_AAseq(protein_filepath)
#     num_aa = 3 # this is the num of amino acids we look ahead in our frame

#     for g in fwd_guides: 
#         # mutates all residues according to the mode, every combination of residue mutations
#         original = g[0][:12] # a string
#         guide_window = g[0][window[0]-1:window[1]] # a string
#         mutateds = [g[0][:window[0]-1] + m + g[0][window[1]:12] for m in make_mutations(guide_window, edit_from, edit_to)] # list of strings 

#         # compares the residues to find which amino acids were altered and catalogs them
#         edits, edit_inds = [], [] # lists of lists, of all edits for all possible mutations
#         start = (-1*g[1])+3
#         orig = original[start:start+(num_aa*3)]
#         for m in mutateds: 
#             # for each possible mutation, come up with the list of amino acid changes
#             edit, edit_ind = find_aa_edits_fwd(m, g, start, orig, num_aa, amino_acid_seq)
#             edits.append(edit)
#             edit_inds.append(edit_ind)
#         # append all information to dataframe
#         g.extend([edits, edit_inds, 'fwd'])
#         assert(len(g)) == 7
        
#     for g in rev_guides: 
#         # mutates all residues according to the mode, every combination of residue mutations
#         original = g[0][:12] # a string
#         guide_window = g[0][window[0]-1:window[1]] # a string
#         mutateds = [g[0][:window[0]-1] + m + g[0][window[1]:12] for m in make_mutations(guide_window, edit_from, edit_to)] # a list of strings

#         # compares the residues to find which amino acids were altered and catalogs them
#         edits, edit_inds = [], [] # lists of lists, of all edits for all possible mutations
#         start = g[1]+1
#         orig = rev_complement(complements, original[start:start+(num_aa*3)])
#         for m in mutateds: 
#             # for each possible mutation, come up with the list of amino acid changes
#             edit, edit_ind = find_aa_edits_rev(m, g, start, orig, num_aa, amino_acid_seq)
#             edits.append(edit)
#             edit_inds.append(edit_ind)
#         # append all information to dataframe
#         g.extend([edits, edit_inds, 'rev'])
#         assert(len(g)) == 7
        
# #     print(pd.DataFrame(fwd_guides + rev_guides))
#     return pd.DataFrame(fwd_guides + rev_guides)

# # function to find all aa edit in the fwd direction
# # given 1 mutant and information on where to start translating
# def find_aa_edits_fwd(m, g, start, orig, num_aa, amino_acid_seq): 
#     mutated = m[start:start+(num_aa*3)] # the mutated string that we are converting
#     original_aa, mutated_aa = '', ''
#     edit, edit_ind = [], []
#     # for the next num_aa amino acids, log each translation and log the changes if there are any
#     for i in range(num_aa): 
#         if not orig[i*3:(i+1)*3].isupper(): # make sure bps are all coding (caps)
#             continue
#         # add translations to a string
#         original_aa += DNA_AA_map[orig[i*3:(i+1)*3]]
#         mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]
#         assert amino_acid_seq[int((g[2]+1+start+(i*3))/3)+1]==original_aa[-1] # check we referenced correct aa
#         # log all changed aa
#         if original_aa[-1] != mutated_aa[-1]: 
#             edit.append(original_aa[-1] + ">" + mutated_aa[-1])
#             edit_ind.append(int((g[2]+1+start+(i*3))/3)+1) # +1 for indexing starting at 1
            
#     if len(edit) == 0: 
#         edit.append('No Change')
#     return edit, edit_ind
    
# # function to find all aa edit in the fwd direction
# # given 1 mutant and information on where to start translating
# def find_aa_edits_rev(m, g, start, orig, num_aa, amino_acid_seq): 
#     mutated = rev_complement(complements, m[start:start+(num_aa*3)])
#     original_aa, mutated_aa = '', ''
#     edit, edit_ind = [], []
#     # for the next num_aa amino acids, log each translation and log the changes if there are any
#     for i in range(num_aa): 
#         if not orig[i*3:(i+1)*3].isupper(): # make sure bps are all coding (caps)
#             continue
#         # add translations to a string
#         original_aa += DNA_AA_map[orig[i*3:(i+1)*3]]
#         mutated_aa += DNA_AA_map[mutated[i*3:(i+1)*3]]
#         assert amino_acid_seq[int((g[2]-10+1+(i*3))/3)+1]==original_aa[-1] # check we referenced correct aa
#         # log all changed aa
#         if original_aa[-1] != mutated_aa[-1]: 
#             edit.append(original_aa[-1] + ">" + mutated_aa[-1])
#             edit_ind.append(int((g[2]-10+1+(i*3))/3))
            
#     if len(edit) == 0: 
#         edit.append('No Change')
#     return edit, edit_ind

def annotate_df(input_df, exonDicts, exonPosDicts, cdsDict, domainDict_all, lenDict):
    # Read in guides, genomic positions, other required annotations.
    # Store in pandas data frame. Each annotation produced by this script will be
    # added as a column to annotated_df.
    annotated_df = input_df.copy()

    #%% Annotate genomic position of editing window for each guide (INCLUSIVE)

    def annot_window(guidePos, guideSense, geneStrand):
        
        # Return (-1,-1) for controls
        if guidePos < 0:
            return (-1,-1)
        
        # If gene on + strand, sense = fwd, antisense = rev
        elif geneStrand == 'plus':
            if guideSense == 'sense':
                return (guidePos+3,guidePos+7)
            elif guideSense == 'antisense':
                return (guidePos-7,guidePos-3)
        
        # If gene on - strand, sense = rev, antisense = fwd
        elif geneStrand == 'minus':
            if guideSense == 'sense':
                return (guidePos-7,guidePos-3)
            elif guideSense == 'antisense':
                return (guidePos+3,guidePos+7)

    annotated_df = annotated_df.assign(Edit_window=annotated_df.apply(
        lambda x: annot_window(x.sgRNA_pos, x.sgRNA_strand, x.Gene_strand), axis=1))

    #%% Annotate the targeted exon (based on first position of protospacer)

    def annot_exon(guidePos, geneID, exonPosDicts):
        
        # Return -1 for controls
        if guidePos < 0:
            return -1
        else:
            # Iterate over each exon (+flanking introns) and check if guide is inside
            for exon, exonPos in exonPosDicts[geneID].items():
                if (guidePos >= exonPos[0]) and (guidePos <= exonPos[1]):
                    return exon
                elif (exon == len(exonPosDicts[geneID])-1):
                    print('Encountered a guide that does not fall within any exon.')
                    return 'None'
        
    annotated_df = annotated_df.assign(Targeted_exon=annotated_df.apply(
        lambda x: annot_exon(x.sgRNA_pos, x.Gene, exonPosDicts), axis=1))

    #%% Annotate window overlap
    # Label whether guides (besides controls) target exonic or intronic regions.
    # Exon if editing window completely within exon. If completely outside exon,
    # label as Intron, 5'-UTR, or 3'-UTR as appropriate. If window spans exon-
    # intron junction, label as Exon/Intron, Exon/5'-UTR, or Exon/3'-UTR.

    def annot_overlap(exon_num, window, geneStrand, geneID, exonPosDicts):
            
        # Annotate controls
        if exon_num < 0:
            return 'Control'
        
        # Make a list of editing window genomic positions
        window_list = np.arange(window[0], window[1]+1, 1)

        # Get genomic coordinates of exon (inclusive)
        # Note that if gene is on minus strand, exon_i is actually end of exon
        exon_i = exonPosDicts[geneID][exon_num][0] + pad
        exon_f = exonPosDicts[geneID][exon_num][1] - pad
        # Make a list of all genomic positions within exon
        exon_list = np.arange(exon_i, exon_f+1, 1)
        
        # Check if there is complete overlap between editing window and exon
        if all(pos in exon_list for pos in window_list):
            return 'Exon'
        # Check if there is partial overlap between editing window and exon
        elif any(pos in exon_list for pos in window_list):
            if exon_num == 0:
                if (geneStrand == 'plus') and (window[0] < exon_i):
                    return "Exon/5'-UTR"
                elif (geneStrand == 'minus') and (window[1] > exon_f):
                    return "Exon/5'-UTR"
                else:
                    return 'Exon/Intron'
            elif exon_num == len(exonPosDicts[geneID])-1:
                if (geneStrand == 'plus') and (window[1] > exon_f):
                    return "Exon/3'-UTR"
                elif (geneStrand == 'minus') and (window[0] < exon_i):
                    return "Exon/3'-UTR"
                else:
                    return 'Exon/Intron'
            else:
                return 'Exon/Intron'
        # If neither above condition is true, window is completely outside exon
        else:
            if exon_num == 0:
                if (geneStrand == 'plus') and (window[0] < exon_i):
                    return "5'-UTR"
                elif (geneStrand == 'minus') and (window[0] > exon_f):
                    return "5'-UTR"
                else:
                    return 'Intron'
            elif exon_num == len(exonPosDicts[geneID])-1:
                if (geneStrand == 'plus') and (window[0] > exon_f):
                    return "3'-UTR"
                elif (geneStrand == 'minus') and (window[0] < exon_i):
                    return "3'-UTR"
                else:
                    return 'Intron'
            else:
                return 'Intron'

    annotated_df = annotated_df.assign(Win_overlap=annotated_df.apply(
        lambda x: annot_overlap(x.Targeted_exon, x.Edit_window,
                                x.Gene_strand, x.Gene, exonPosDicts), axis=1))

    #%% Annotate with number of 'C's in editing window and whether C is present

    annotated_df = annotated_df.assign(C_count=annotated_df.apply(
        lambda x: x.sgRNA_seq[3:8].count('C'), axis=1))

    def annot_isC(C_count):
        
        #Label whether there is at least one C in editing window
        if C_count >= 1:
            return 'C'
        elif C_count == 0:
            return 'No_C'
        else: print('This should not occur')

    annotated_df = annotated_df.assign(is_C=annotated_df.apply(
        lambda x: annot_isC(x.C_count), axis=1))

    #%% Annotate splice site targeting guides
    # Four categories: 'Splice-donor', 'Splice-acceptor', 'None', and 'Control'
    # 'Splice-donor' means guide edits the 'GT' at 5' of intron
    # 'Splice-acceptor' means guide edits the 'AG' at 3' of intron
    # For C editor, only antisense guides can qualify. Strategy here is to check
    # (1) that guide is antisense, (2) that its editing window overlaps with
    # the 'C' at either splice site within that exon, (3) that there is a 'C' there.

    def annot_splice(winOverlap, guideStrand, geneStrand, window, geneID, exon,
                    exonDicts, exonPosDicts):
        
        # Annotate controls
        if winOverlap == 'Control':
            return 'Control'
        else:
            if guideStrand == 'antisense':
                # Define exon boundaries (stored as tuple)
                exon_coords = exonPosDicts[geneID][exon]
                # Identify genomic location of splice site editable bases
                if geneStrand == 'plus':
                    donor_loc = exon_coords[1]-(pad-1)
                    acceptor_loc = exon_coords[0]+(pad-1)
                elif geneStrand == 'minus':
                    donor_loc = exon_coords[0]+(pad-1)
                    acceptor_loc = exon_coords[1]-(pad-1)
                
                # Check that there is an editable base there
                # Define exon+flanking intron sequence (given 5' to 3' along gene)
                exon_seq = str(exonDicts[geneID][exon]).upper()
                donor_present = ('G' == exon_seq[-pad])
                acceptor_present = ('G' == exon_seq[(pad-1)])
                
                # Expand editing window boundaries into a list of all positions
                window_list = np.arange(window[0], window[1]+1, 1)
                
                # Initialize string to hold annotation (doing it this way just in
                # case somehow a guide hits both an acceptor and donor for small exon)
                splice_note = ''
                
                # Check for donor editing
                if all([donor_loc in window_list, donor_present,
                        exon!=(len(exonPosDicts[geneID])-1)]):
                    splice_note += 'Splice-donor'
                # Check for acceptor editing
                if all([acceptor_loc in window_list, acceptor_present, exon!=0]):
                    if len(splice_note) == 0:
                        splice_note += 'Splice-acceptor'
                    else:
                        splice_note += '/Splice-acceptor'
                
                # Return annotation
                if len(splice_note) == 0:
                    return 'None'
                else:
                    return splice_note
            else:
                return 'None'

    annotated_df = annotated_df.assign(Splice_check=annotated_df.apply(
        lambda x: annot_splice(x.Win_overlap, x.sgRNA_strand, x.Gene_strand,
                            x.Edit_window, x.Gene, x.Targeted_exon,
                            exonDicts, exonPosDicts), axis=1))

    #%% Annotate residue(s) and CDS position(s) targeted by each guide

    # Process each list of exon lengths in lenDict into a list of CDS cumulative
    # nucleotide numbers (stored in dictionary with key being gene as before)
    shiftDict = {}
    for gene, lenList in lenDict.items():
        shiftList = []
        for i in range(len(lenList)):
            shiftList.append(sum(lenList[:i]))
        shiftDict[gene] = shiftList
        
    def annot_target_CDS(winOverlap, window, exon_num, geneStrand, geneID,
                        exonPosDicts, shiftDict):
        
        # Assign (-1,-1) to all guides whose windows do not overlap exons
        if not 'Exon' in winOverlap:
            return (-1,-1)
        
        #i and f refer to lower and upper boundaries in genomic coordinates
        #Note if gene is along minus strand, i is downstream of f along the gene
        window_i = window[0]
        window_f = window[1]
        exon_i = exonPosDicts[geneID][exon_num][0] + pad
        exon_f = exonPosDicts[geneID][exon_num][1] - pad
        if window_f > exon_f:
            window_f = exon_f
        if window_i < exon_i:
            window_i = exon_i
        
        #Find targeted interval along CDS (1-indexed, inclusive boundaries)
        if geneStrand == 'plus':
            CDS_i = shiftDict[geneID][exon_num] + (window_i - exon_i + 1)
            CDS_f = shiftDict[geneID][exon_num] + (window_f - exon_i + 1)
        elif geneStrand == 'minus':
            CDS_i = shiftDict[geneID][exon_num] + (exon_f - window_f + 1)
            CDS_f = shiftDict[geneID][exon_num] + (exon_f - window_i + 1)

        return (CDS_i, CDS_f)

    annotated_df = annotated_df.assign(Target_CDS=annotated_df.apply(
        lambda x: annot_target_CDS(x.Win_overlap, x.Edit_window, x.Targeted_exon,
                                x.Gene_strand, x.Gene, exonPosDicts,
                                shiftDict), axis=1))

    def annot_target_resi(targetCDS):
        
        # Only consider exonic guides
        if targetCDS == (-1,-1):
            return -1
        
        CDS_i = targetCDS[0]
        CDS_f = targetCDS[1]
        
        aa_i = math.ceil(CDS_i / 3)
        aa_f = math.ceil(CDS_f / 3)
        
        return (aa_i, aa_f)

    annotated_df = annotated_df.assign(Target_resi=annotated_df.apply(
        lambda x: annot_target_resi(x.Target_CDS), axis=1))

    #%% Annotate mutation type for exonic guides with C

    def annot_mutations(guideID, winOverlap, num_C, window, guideSense, geneID,
                        spliceCheck, cdsDict):
        # Returns tuples (guideID, Mut_type, Mut_list_all, Mut_list)
        
        # Annotate controls
        if winOverlap == 'Control':
            return (guideID, 'Control', 'None', 'None')
        
        # Annotate guides without a C in the editing window
        if num_C == 0:
            if 'Exon' in winOverlap:
                return (guideID, 'No_C/Exon', 'None', 'None')
            else:
                return (guideID, 'No_C/Non-exon', 'None', 'None')
            
        # Annotate strictly non-exonic guides (not overlapping exons at all)
        if not 'Exon' in winOverlap:
            if 'Splice' in spliceCheck:
                return (guideID, 'Splice', 'None', 'None')
            else:
                return (guideID, 'Non-exon', 'None', 'None')
            
        # Window is the editing window boundaries, inclusive, 1-indexed
        window_seq = cdsDict[geneID][(window[0]-1):(window[1])]
        
        if guideSense == 'sense':
            base = 'C'
            new_base = 'T'
        elif guideSense == 'antisense':
            base = 'G'
            new_base = 'A'
        
        # All remaining guides (1) contain a C in the editing window and
        # (2) are either strictly Exonic or span an Exon/Non-exon junction
        
        # Thus, if there isn't an editable base within the exonic portion of the
        # window, the guide must span a junction, and the editable base must
        # be in the intronic part, so we annotate as Non-Exon (unless splice).
        if base not in window_seq:
            if 'Splice' in spliceCheck:
                return (guideID, 'Splice', 'None', 'None')
            else:
                return (guideID, 'Non-exon', 'None', 'None')
        
        # Generate list of all possible mutated versions of window_seq
        mutList = generateMutCom(window_seq,
                                generateSubsets(find_Indexes(base, window_seq, [])),
                                new_base)
        
        # Identify fully mutated window_seq (all C>T or G>A within window)
        maxMut = generateMutCom(window_seq,
                                [tuple(find_Indexes(base, window_seq, []))],
                                new_base)
        maxMut = maxMut[0] #Unpack from list
        
        #Assemble the codons surrounding the editing site
        window_0 = (window[0]-1, window[1]-1) #Convert editing window to zero-index
        first_codon_i = window_0[0] - window_0[0]%3
        last_codon_i = (window_0[1] + 3) - window_0[1]%3
        first_codon_adj = cdsDict[geneID][first_codon_i:window_0[0]]
        last_codon_adj = cdsDict[geneID][window_0[1]+1:last_codon_i]
        CDS_codonstr = first_codon_adj + window_seq + last_codon_adj #Done
        
        #Assemble analogous codons for potential mutated copies
        mut_codonstr = []
        for mut in mutList:
            mut_codonstr.append(first_codon_adj + mut + last_codon_adj)
        maxMut_str = first_codon_adj + maxMut + last_codon_adj #fully mutated
        
        #Break apart into individual codons
        CDS_codons = [CDS_codonstr[i:i+3] for i in range(0, len(CDS_codonstr), 3)]
        
        #Iterate through list of all possible combinations of base changes and
        #compile list of non-synonymous mutations
        all_mutlist = []
        for mut_seq in mut_codonstr:
            mut_codons = [mut_seq[i:i+3] for i in range(0, len(mut_seq), 3)]
            
            singlecom_mutlist = []
            for i, codon in enumerate(CDS_codons):
                aa_index = int((window_0[0]-window_0[0]%3)/3 + 1 + i) #1-index
                aa_CDS = get_aa(codon)
                aa_mut = get_aa(mut_codons[i])
                if aa_CDS == aa_mut:
                    continue
                else:
                    mut_note = str(aa_CDS+str(aa_index)+aa_mut)
                    singlecom_mutlist.append(mut_note)
            
            if len(singlecom_mutlist) > 0:
                all_mutlist.extend(singlecom_mutlist)
        
        #Make analogous list for fully mutated
        max_all_mutlist = []
        max_mut_codons = [maxMut_str[i:i+3] for i in range(0, len(maxMut_str), 3)]
        for i, codon in enumerate(CDS_codons):
            aa_index = int((window_0[0]-window_0[0]%3)/3 + 1 + i) #1-index
            aa_CDS = get_aa(codon)
            aa_mut = get_aa(max_mut_codons[i])
            if aa_CDS == aa_mut:
                continue
            else:
                mut_note = str(aa_CDS+str(aa_index)+aa_mut)
                max_all_mutlist.append(mut_note)
        
        #Prep annotation for all possible mutations (combinations of base edits)
        if len(all_mutlist) > 0:
            final_all_mutlist = set(all_mutlist)
            all_annot = ''
            for index, mut in enumerate(final_all_mutlist):
                if index == 0:
                    temp = mut
                else: temp = ', ' + mut
                all_annot += temp
        else:
            all_annot = 'None'
        
        # Prep fully edited annotation, including mutation type, and produce output
        if len(max_all_mutlist) > 0:
            final_list = ''
            nonsense_check = False
            for index, mut in enumerate(max_all_mutlist):
                if '*' in mut:
                    nonsense_check = True
                if index == 0:
                    temp = mut
                else: temp = ', ' + mut
                final_list += temp
            if nonsense_check:
                return (guideID, 'Nonsense', all_annot, final_list)
            else:
                if 'Splice' in spliceCheck:
                    return (guideID, 'Splice', all_annot, final_list)
                else:
                    return (guideID, 'Missense', all_annot, final_list)
        else:
            if 'Splice' in spliceCheck:
                return (guideID, 'Splice', all_annot, 'None')
            else:
                return (guideID, 'Silent', all_annot, 'None')

    mut_notes = annotated_df.apply(lambda x: annot_mutations(x.sgRNA_ID, x.Win_overlap, x.C_count,
                                x.Target_CDS, x.sgRNA_strand, x.Gene,
                                x.Splice_check, cdsDict), axis=1)
    mut_notes_df = pd.DataFrame(mut_notes.tolist(),
                                columns=['sgRNA_ID','Mut_type','Mut_list_all','Mut_list'])
    annotated_df = annotated_df.merge(mut_notes_df, on='sgRNA_ID')

    #%% Assign each guide a position (residue numbering)
    # All guides apart from exonic ones: set to -1
    # Exonic guides: assign position based on 1-indexed protein sequence
        # No C: Use the center of the editing window (convert to residue number)
        # Yes C:
            # Silent: Use the center of the editing window
            # Missense: Take mutated residue numbers and average
            # Nonsense: Take nonsense mutation residue numbers (avg if multiple)

    def annot_editsite(window, mutType, mutList):
        
        # Set edit site to -1 if non-exonic or control
        if mutType in ['Non-exon', 'No_C/Non-exon', 'Control', 'Splice']:
            return -1
        elif (mutType == 'No_C/Exon') or (mutType == 'Silent'):
            avg_site = sum(window) / len(window)
            # If window contains odd number of positions, avg site is discrete base
            if (window[1]-window[0])%2 == 0:
                return math.ceil(avg_site / 3)
            else:
                # Otherwise avg site is between two bases. Find residues each maps to
                upper_site = math.ceil(avg_site)
                lower_site = math.floor(avg_site)
                upper_aa = math.ceil(upper_site/3)
                lower_aa = math.ceil(lower_site/3)
                # If both bases map to same residue, assign to that residue
                if upper_aa == lower_aa:
                    return upper_aa
                # Otherwise take average
                else:
                    return (upper_aa+lower_aa)/2
        elif (mutType == 'Nonsense'):
            muts_parsed = mutList.split(', ')
            resi_nums = [int(mut[1:-1]) for mut in muts_parsed if '*' in mut]
            return sum(resi_nums)/len(resi_nums)
        elif (mutType == 'Missense'):
            muts_parsed = mutList.split(', ')
            resi_nums = [int(mut[1:-1]) for mut in muts_parsed]
            resi_nums_unique = set(resi_nums)
            return sum(resi_nums_unique)/len(resi_nums_unique)
        else:
            print('This should not occur!')

    annotated_df = annotated_df.assign(Edit_site=annotated_df.apply(
        lambda x: annot_editsite(x.Target_CDS, x.Mut_type, x.Mut_list), axis=1))

    #%% Annotate domain

    # Label all exonic guides by domain (e.g. any guides with an Edit_site not -1).
    # For non-exonic guides, label as Control, Splice, Intron, 5'-UTR, or 3'-UTR.
    # Exonic guides that do not fall within a provided domain are labeled 'Undefined'.
    def annot_domain(site, winOverlap, mutType, geneID, domainDict_all):
        
        # Annotate controls
        if mutType == 'Control':
            return 'Control'
        # Annotate splice-site targeting guides
        elif mutType == 'Splice':
            return 'Splice'
        # Annotate non-exonic guides
        elif 'Non-exon' in mutType:
            # Assign Intron, 5'-UTR, or 3'-UTR as appropriate
            if "5'-UTR" in winOverlap:
                return "5'-UTR"
            elif "3'-UTR" in winOverlap:
                return "3'-UTR"
            elif 'Intron' in winOverlap:
                return 'Intron'
            else:
                print('This should not occur!')
                return
        # If none of above conditions is true, the guide must be exonic
        else:
            # Use gene annotation to pull out correct domain boundary dictionary
            bounds_Dict = domainDict_all[geneID]
            # Iterate over domains and check if guide lies within each one
            for domain, bounds in bounds_Dict.items():
                if bounds[0] <= site <= bounds[1]:
                    return domain
                else:
                    continue
            # If guide was not found in any domain, label as 'Undefined'
            return 'Undefined'

    annotated_df = annotated_df.assign(Domain=annotated_df.apply(
        lambda x: annot_domain(x.Edit_site, x.Win_overlap, x.Mut_type, x.Gene,
                            domainDict_all), axis=1))

    #%% Annotate plot sites
    # Assign a pseudocount for plotting to any guides that have Edit_site = -1.
    # These are guides with Mut_type in {Non-exon, No_C/Non-exon, Splice, Control}.
    # For guides that have been assigned a true Edit_site (guides with Mut_type in
    # {No_C/Exon, Silent, Missense, Nonsense}, this will copy that annotation.

    # Set the spacings between pseudocounts. Resulting scatter plot x-axis:
    # Exonic // Non-exonic // Splice // Controls: Intergenic / Non-targeting / Essential) 
    pseudospacer = 25 # Sets // spacing
    ctrl_spacer = 10 # Sets / spacing
    pseudodensity = 0.5 # Sets spacing between adjacent pseudocounts
    total_aas = 472 # Max number that Edit_site can take on
    def annot_plotsite(df_input, pseudospacer, ctrl_spacer, pseudodensity, total_aas):
        
        # Count how many guides are in each category
        counts = df_input['Mut_type'].value_counts()
        if ('Non-exon' in counts.index.values) & ('No_C/Non-exon' in counts.index.values):
            n_nonexon = counts['Non-exon']+counts['No_C/Non-exon']
        else: n_nonexon = 0
        
        if 'Splice' in counts.index.values:
            n_splice = counts['Splice']
        else: n_splice = 0
            
        if 'Control' in counts.index.values:
            n_control = counts['Control']
        else: n_control = 0
        
        # Control guide subcounts
        ctrl_counts = df_input['Gene'].value_counts()
        if 'NON-GENE' in ctrl_counts.index.values:
            n_nongene = ctrl_counts['NON-GENE']
        else: n_nongene = 0
        if 'NO_SITE' in ctrl_counts.index.values:
            n_nongene = ctrl_counts['NO_SITE']
        else: n_nosite = 0
        if 'ESSENTIAL' in ctrl_counts.index.values:
            n_essential = ctrl_counts['ESSENTIAL']
        else: n_essential = 0
        
        # Checks
        if n_control != n_nongene+n_nosite+n_essential:
            print('Control guides do not add up!')
        
        # Create dataframe to store values
        df_store = df_input[['sgRNA_ID','Edit_site']].copy()
        df_store = df_store.set_index('sgRNA_ID')
        df_store['Plot_site'] = 0
        
        # Initialize
        i_nonexon = total_aas + pseudospacer
        i_splice = i_nonexon + n_nonexon*pseudodensity + pseudospacer
        i_nongene = i_splice + n_splice*pseudodensity + pseudospacer
        i_nosite = i_nongene + n_nongene*pseudodensity + ctrl_spacer
        i_essential = i_nosite + n_nosite*pseudodensity + ctrl_spacer
        
        for row in df_input.itertuples():
            
            guideID = row.sgRNA_ID
            
            if row.Edit_site < 0:
                if (row.Mut_type == 'Non-exon') or (row.Mut_type == 'No_C/Non-exon'):
                    df_store.loc[guideID, 'Plot_site'] = i_nonexon
                    i_nonexon += pseudodensity
                elif row.Mut_type == 'Splice':
                    df_store.loc[guideID, 'Plot_site'] = i_splice
                    i_splice += pseudodensity
                elif row.Mut_type == 'Control':
                    if row.Gene == 'NON-GENE':
                        df_store.loc[guideID, 'Plot_site'] = i_nongene
                        i_nongene += pseudodensity
                    elif row.Gene == 'NO_SITE':
                        df_store.loc[guideID, 'Plot_site'] = i_nosite
                        i_nosite += pseudodensity
                    elif row.Gene == 'ESSENTIAL':
                        df_store.loc[guideID, 'Plot_site'] = i_essential
                        i_essential += pseudodensity
                    else:
                        print('This should not occur!')
                else:
                    print('This should not occur!')
            else:
                df_store.loc[guideID, 'Plot_site'] = row.Edit_site
            
        df_store = df_store.reset_index()
        df_store = df_store.drop(columns='Edit_site')
        
        return df_store
        

    # NOTE: Only annotate plot sites if there's only one gene
    if plot_site_check:
        df_plotsites = annot_plotsite(annotated_df, pseudospacer, ctrl_spacer,
                                    pseudodensity, total_aas)
        annotated_df = annotated_df.merge(df_plotsites, on='sgRNA_ID')

    annotated_df.to_csv(out_prefix+'_Annotated_Library_Full_AR.csv', index=False)

    #%% Base editing heatmap
    # Plot base editing scope heatmap.

    list_aas = ['R','H','K','D','E',
                'S','T','N','Q','C','G','P',
                'A','V','I','L','M','F','Y','W','*']

    df_muts = pd.DataFrame(np.zeros((21,21)), index=list_aas, columns=list_aas)

    # Make a list of all possible codons or trinucleotides
    list_codons = list(itertools.product('ATCG', repeat=3))
    list_codons = [''.join(tup) for tup in list_codons]

    for codon in list_codons:
        
        # AA that the codon encodes
        start_aa = get_aa(codon)
        
        # Generate all possible C to T mutations
        C_subsets = generateSubsets(find_Indexes('C', codon, []))
        C_muts = generateMutCom(codon, C_subsets, 'T')
        
        # Generate all possible G to A mutations
        G_subsets = generateSubsets(find_Indexes('G', codon, []))
        G_muts = generateMutCom(codon, G_subsets, 'A')
        
        # All possible mutations
        all_muts = C_muts + G_muts
        
        # All possible outcome codons
        end_aas = [get_aa(mut_codon) for mut_codon in all_muts]
        
        # Fill in df_muts
        for result in end_aas:
            df_muts.loc[start_aa, result] = 0.3

    # Set silent mutation squares to zero
    for aa in list_aas:
        df_muts.loc[aa,aa] = 0    

    def plot_heatmap(df_plot, title, out_prefix):
        
        # Set plotting parameters
        sns.set_context('talk')
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        sns.set(font_scale=1.2)
        sns.set_style('whitegrid')
        
        #cmap = sns.diverging_palette(240, 10, l=10, as_cmap=True)
        
        # Generate heatmap
        fig, ax = plt.subplots(figsize=(10,10))
        sns.heatmap(data=df_plot, ax=ax, square=True, cmap='Reds', vmin=0, vmax=1,
                    linewidths=0.5, linecolor='gray',
                    xticklabels=True, yticklabels=True, cbar_kws={"shrink": .70})
        for edge,spine in ax.spines.items():
            spine.set_visible(True)
            spine.set_color('k')
            spine.set_linewidth(2)
        plt.xticks(rotation=0, horizontalalignment='center')
        plt.yticks(rotation=0, horizontalalignment='center')
        plt.xlabel('Mutant residue')
        plt.ylabel('Initial residue')
        plt.title(title)
        line_locs = [3,5,9,12,17,20]
        for loc in line_locs:
            ax.axhline(y=loc, color='k', lw=2)
            ax.axvline(x=loc, color='k', lw=2)
        
        plt.savefig(out_prefix+'_heatmap.pdf', format='pdf', bbox_inches='tight')
        plt.close()

    # Make heatmap
    plot_heatmap(df_muts, 'Possible Non-Silent Mutations with C to T Base Editor',
                out_prefix+'_BE_mutation')


    # Plot gene-specific heatmap

    def library_heatmap(df_annot, df_plot):
        
        df_out = df_plot.copy()
        
        # Isolate mutation lists for missense and nonsense guides and combine
        df_missense = df_annot.loc[df_annot['Mut_type']=='Missense'].copy()
        df_nonsense = df_annot.loc[df_annot['Mut_type']=='Nonsense'].copy()
        missense_list = df_missense['Mut_list'].tolist()
        nonsense_list = df_nonsense['Mut_list'].tolist()
        total_list = missense_list + nonsense_list
        
        # Parse entries in mutation lists and set values in dataframe for heatmap
        for mutString in total_list:
            mutList = mutString.split(', ')
            for mut in mutList:
                aa_i = mut[0]
                aa_f = mut[-1]
                df_out.loc[aa_i,aa_f] = 0.8
        
        return df_out
        
    df_librarymuts = library_heatmap(annotated_df, df_muts)

    # Make heatmap
    plot_heatmap(df_librarymuts, 'AR Library Mutations, C to T Base Editor',
                out_prefix+'_AR_library')

    #%% Miscellaneous plots

    # Barplot by Mut_type

    def plot_mutTypes(df_annot, title, out_prefix):
        
        df_data = df_annot.copy()
        
        # Set order that mutation types will appear
        type_list = ['Nonsense','Missense','Silent','No_C/Exon','Non-exon',
                    'No_C/Non-exon','Splice','Control']
        
        # ColorBrewer2, 9 data classes, qualitative, 4th color scheme hex codes
        color_list = ['#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', '#d9d9d9',
                    '#8dd3c7', '#ffffb3', '#bebada']
        
        # Set plotting parameters
        sns.set_context('talk')
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        sns.set_palette('bright')
        sns.set_style('white')
        sns.set_style('ticks')
        
        # Make plot
        fig, ax = plt.subplots(figsize=(10,6))
        g = sns.countplot(y='Mut_type', data=df_data, order=type_list, hue_order=type_list,
                        palette=color_list, ax=ax, saturation=1, edgecolor='black',
                        linewidth=1.5)
        for p in g.patches:
            guide_count = p.get_width()
            g.text(guide_count+1, p.get_y()+p.get_height()/2, guide_count,
                ha='left', va='center')
        #plt.tight_layout()
        plt.title(title)
        sns.despine()
        
        # Save
        plt.savefig(out_prefix+'_Mut_type_barplot.pdf', format='pdf', bbox_inches='tight')
        plt.close()

    plot_mutTypes(annotated_df, 'AR Library', out_prefix+'_AR_library')

    # Histogram of distances between guides

    def plot_histogram(df_annot, title, out_prefix):
        
        df_data = df_annot.copy()
        
        # Retain guides with true Edit sites (i.e. not -1)
        df_data = df_data.loc[df_data['Mut_type']=='Missense']
        df_data = df_data.sort_values(by='Edit_site')
        list_sites = df_data['Edit_site'].tolist()
        
        # Find distances between each guide
        list_dists = []
        for i in range(1, len(list_sites)):
            list_dists.append(list_sites[i]-list_sites[i-1])
        
        # Set plotting parameters
        sns.set_context('talk')
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
        sns.set_palette('deep')
        sns.set_style('ticks')
        
        # Plot
        sns.distplot(list_dists, kde=False, color='#80b1d3',
                    hist_kws={'alpha':1,'edgecolor':'black'})
        sns.despine()
        plt.tight_layout()
        plt.xlabel('Distance between guides')
        plt.ylabel('Frequency')
        plt.title(title)
        
        # Save
        plt.savefig(out_prefix+'_missense_distance_histogram.pdf', format='pdf',
                    bbox_inches='tight')
        plt.close()

    plot_histogram(annotated_df, 'AR Missense guides', out_prefix+'_AR')

    #%% Identify guides whose protospacers+PAMs are within exons or span boundaries.
    # Unlike above annotations, based on whole protospacer+PAM not just edit window
    # Need this information so we can design DNMT3A CDS that is not recognized by
    # these guides.

    def annot_exonoverlap(guideID, guidePos, guideSense, geneStrand, exon_num, geneID, exonPosDicts):
        
        # Return No for controls
        if guidePos < 0:
            return (guideID, 'No', 0)
        
        # Now figure out the coordinates of the protospacer+PAM
        
        # If gene on + strand, sense = fwd, antisense = rev
        if geneStrand == 'plus':
            if guideSense == 'sense':
                coords = (guidePos, guidePos+22)
            elif guideSense == 'antisense':
                coords = (guidePos-22, guidePos)
        
        # If gene on - strand, sense = rev, antisense = fwd
        elif geneStrand == 'minus':
            if guideSense == 'sense':
                coords = (guidePos-22, guidePos)
            elif guideSense == 'antisense':
                coords = (guidePos, guidePos+22)
        
        # Make a list of all guide positions
        guide_list = np.arange(coords[0], coords[1]+1, 1)
        
        # Retrieve exon CDS boundaries
        # Get genomic coordinates of exon (inclusive)
        # Note that if gene is on minus strand, exon_i is actually end of exon
        exon_i = exonPosDicts[geneID][exon_num][0] + pad
        exon_f = exonPosDicts[geneID][exon_num][1] - pad
        # Make a list of all genomic positions within exon
        exon_list = np.arange(exon_i, exon_f+1, 1)
        
        pos_sum = sum(pos in exon_list for pos in guide_list)
        if all(pos in exon_list for pos in guide_list):
            return (guideID, 'Exon', pos_sum)
        elif any(pos in exon_list for pos in guide_list):
            return (guideID, 'Boundary', pos_sum)
        else:
            return (guideID, 'No', pos_sum)
        

    ex_overlap_notes = annotated_df.apply(
        lambda x: annot_exonoverlap(x.sgRNA_ID, x.sgRNA_pos, x.sgRNA_strand, x.Gene_strand,
                                    x.Targeted_exon, x.Gene, exonPosDicts), axis=1)
    ex_overlap_notes_df = pd.DataFrame(ex_overlap_notes.tolist(),
                                    columns=['sgRNA_ID','Exon_overlap','Exon_overlap_num'])
    annotated_df = annotated_df.merge(ex_overlap_notes_df, on='sgRNA_ID')

    return annotated_df

def export_csvs(annotated_df, out_prefix):
    # Export csvs
    within_exon_df = annotated_df.loc[annotated_df['Exon_overlap']=='Exon'].copy()
    boundary_df = annotated_df.loc[annotated_df['Exon_overlap']=='Boundary'].copy()
    non_exon_df = annotated_df.loc[annotated_df['Exon_overlap']=='No'].copy()

    key_cols = ['sgRNA_ID', 'sgRNA_seq', 'sgRNA_pos', 'sgRNA_strand','Mut_type',
                'Mut_list','Edit_site','Domain','Exon_overlap','Exon_overlap_num']

    within_exon_df[key_cols].to_csv(out_prefix+'_within_exon_guides.csv', index=False)
    boundary_df[key_cols].to_csv(out_prefix+'_exon_intron_boundary_guides.csv', index=False)
    non_exon_df[key_cols].to_csv(out_prefix+'_not_in_exon_guides.csv', index=False)

if __name__ == "__main__":
    #%% For inputs generated by CRISPOR
    CRISPOR_output = '211116_CCNK_Output.tsv' # tsv output with all of the gRNAs generated by CRISPOR
    input_fasta = '211116_CCNK_Input.fasta' # fasta input generated by UCSC genome browser that was piped into CRISPOR (sequence of gene)
    gene_name = 'CCNK' # gene name
    strand_info = 'plus' # strand information on the gene (whether it is on the plus or minus strand in the genome)
    PAM = 'NG' # PAM used (e.g. NGG, NG)

    guide_input_csv = "guideInput_CCNK.csv" # name of the new file that was converted from tsv output that will be the guide input for the script

    #tuples defining domain boundaries (inclusive, 1-indexed resi #s). Any exonic guides not within any of these boundaries will be labeled 'Undefined'.
    domainDict = {'NTD':(1,554),
                'DBD':(555,623),
                'Hinge':(624,664),
                'LBD':(665,919),
                } 

    # Output file name
    out_dir = 'CCNK_Output'
    out_prefix = os.path.join(out_dir,'211116')

    # Padding info
    pad = 20 # padding around exons, if you used CRISPOR, default is usually 20


    # @params fr l to r: Output.tsv = output generated by CRISPOR, Input.fasta = input that was piped into CRISPOR, GENE = gene name, strand info (plus/minus), NG = PAM used
    guideInput = convert_guideInput(CRISPOR_output, input_fasta, gene_name, strand_info, PAM)
    # @params fr l to r: guideInput.csv = name of output file (to be used for the script)
    guideInput.to_csv(guide_input_csv, sep=',',index=False)

    #%% Key inputs

    # Dictionary with keys being possible values of 'Gene' column (excluding
    # controls) and values being the corresponding fasta files
    geneDict = {gene_name:input_fasta}

    # Input annotations. Requires 6 columns: sgRNA_ID, sgRNA_seq, Gene, sgRNA_pos,
    # sgRNA_strand (sense or antisense), Gene_strand (plus or minus)
    # Should be the same as line 39
    guideInput = guide_input_csv

    # Script currently does not support changing base editor type and assumes
    # +4 to +8 editing window and C>T editor
    #editrange = (3,7) #0-index
    #base = 'C'

    # Define domains (this info will be used to assign guides to domains).
    # First, make dictionaries for each gene, where keys are domains and values are
    # tuples defining domain boundaries (inclusive, 1-indexed resi #s). Any exonic
    # guides not within any of these boundaries will be labeled 'Undefined'.
    domainDict1 = domainDict
    #domainDict2 = {'domain1':(start,end),
    #               'domain2':(start,end)}
    # Second, define master domain dictionary with keys being possible genes
    # (make sure these match the keys in geneDict!) and values being above dict(s).
    domainDict_all = {gene_name:domainDict1}

    # Make Plot_site annotation? This creates pseudocount for non-exonic, control,
    # and splice guides to enable them to be plotted. However, this annotation
    # should only be enabled if the library contains only one gene.
    plot_site_check = True

    #%% Analysis set up

    # Read each gene's fasta file and parse exon sequences and genomic positions.
    # Store within dictionaries where keys are genes and values are corresponding
    # parsed information (also dictionaries).
    exonDicts = {}
    exonPosDicts = {}
    for gene, fastaFile in geneDict.items():
        exonDicts[gene], exonPosDicts[gene] = processIFile(fastaFile)

    # Prepare reference coding sequences for each gene (stored in dictionary like above)
    cdsDict = {}
    for gene in geneDict.keys():
        full_CDS = ''
        for i in range(len(exonDicts[gene])):
            seq = str(exonDicts[gene][i])
            firstIntronNuc = len(seq) - pad # because slicing does not include second index
            full_CDS += seq[pad:firstIntronNuc]
        cdsDict[gene] = full_CDS


    # Prepare reference lists of exon lengths  (stored in dictionary like above)
    lenDict = {}
    for gene in geneDict.keys():
        lenDict[gene] = [(len(exonDicts[gene][i]) - (pad*2)) for i in range(len(exonDicts[gene]))]

    input_df = pd.read_csv(guideInput)
    annotated_df = annotate_df(input_df, exonDicts, exonPosDicts, cdsDict, domainDict_all, lenDict)
    export_csvs(annotated_df, out_prefix)