#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 12:20:34 2020

@author: nlclarinet8, ceejaylee

Description: This .py file contains functions to accompany the
NZL10196-Guide-Annotation.py script.
"""

#%% Import Packages

import os
import pathlib
import time
import sys
import csv
from collections import OrderedDict, Counter
from itertools import combinations

import numpy as np
import pandas as pd
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

#%% Functions

# Process .fasta into dictionaries with exons (+flanking introns), number as key
def processIFile(fastaFile):
    #import fasta file for comparison
    #exonDict is a dictionary with key = exon #, value = sequence
    #positionDict is a dictionary with key = exon #, value = (start,end) genomic position, inclusive
    exonDict = {}
    exonPosDict = {}

    for exon in SeqIO.parse(fastaFile,'fasta'): #SeqIO.parse takes a file handle (or filename) and format name, and returns a SeqRecord iterator
        exon_id = parseName(exon.id) #.id returns name of fasta chunk, parseName is a function that splits the name to only retain its numbering
        exon_seq = exon.seq #.seq returns the sequence of the fasta chunk
        exon_pos = parsePos(exon.description)
        exonDict[exon_id] = exon_seq #construct dictionary entry
        exonPosDict[exon_id] = exon_pos #construct dictionary for positions tuples
    
    return exonDict, exonPosDict

# Function to convert CRISPOR output tsv of gRNAs to a input format for this script
# Ideal input document would be unprocessed output file from CRISPOR (as opposed to filtered output from visualization tool), as the former contains information on genomic coordinates of the gene

def convert_guideInput(gRNA_Output, gene_fasta, Gene, strand_info, PAM):
    CRISPOR_df = pd.read_csv(gRNA_Output, sep='\t')
    exon_tuple = processIFile(gene_fasta)
    guideInput = pd.DataFrame(columns=['sgRNA_ID','sgRNA_seq','Gene','sgRNA_pos','sgRNA_strand','Gene_strand'])
    pop_n = len(PAM)
    CRISPOR_df['targetSeq']=CRISPOR_df['targetSeq'].apply(lambda x: x[:-pop_n])
    CRISPOR_df['exon_no'] = CRISPOR_df['#seqId'].str.split(' ',expand=True)[0].str.split('_',expand=True)[4]
    
    CRISPOR_df['sgRNA_strand'] = CRISPOR_df['guideId'].apply(lambda x: 'antisense' if 'rev' in x else 'sense')
    CRISPOR_df['sgRNA_pos'] = [find_sgRNA_pos(CRISPOR_df['targetSeq'].iloc[i],CRISPOR_df['sgRNA_strand'].iloc[i],strand_info,x,exon_tuple) for i,x in enumerate(CRISPOR_df['exon_no'])]
    if strand_info=='minus':
        CRISPOR_df = CRISPOR_df.sort_values(by=['exon_no','sgRNA_pos'],ascending=[False,False]).reset_index(drop=True)
    elif strand_info=='plus':
        CRISPOR_df = CRISPOR_df.sort_values(by=['exon_no','sgRNA_pos'],ascending=[True,True]).reset_index(drop=True)
    else: print("Error! No strand information available.")
    n_guides = len(CRISPOR_df['targetSeq'])
    n_char = len(str(n_guides))
    CRISPOR_df = CRISPOR_df.assign(n_guide=lambda x: x.index.astype(str).str.zfill(n_char))
    CRISPOR_df['sgRNA_ID'] = [Gene+'_'+str(x).zfill(3) for x in range(1,n_guides+1)]
    guideInput['sgRNA_ID'] = CRISPOR_df['sgRNA_ID']
    guideInput['sgRNA_seq'] = CRISPOR_df['targetSeq']
    guideInput['Gene'] = Gene
    guideInput['sgRNA_pos'] = CRISPOR_df['sgRNA_pos']
    guideInput['sgRNA_strand'] = CRISPOR_df['sgRNA_strand']
    guideInput['Gene_strand'] = strand_info
    return(guideInput)


# HELPER FUNCTIONS
# Calculate sgRNA pos using the exon sequence/chromosome number and sequence of the gRNA
def find_sgRNA_pos(sgRNA_seq, sgRNA_strand, gene_strand, exon_no, processed_fasta_tuple):
    refSeq = processed_fasta_tuple[0][int(exon_no)].upper()
    refSeq_no = processed_fasta_tuple[1][int(exon_no)]
    if sgRNA_strand == 'sense':
        pos_in_exon = refSeq.find(sgRNA_seq)
        if pos_in_exon >= 0:
            if gene_strand =='plus':
                return refSeq_no[0] + pos_in_exon
            elif gene_strand == 'minus':
                return refSeq_no[1] - pos_in_exon
            else: print("Error! Gene strand info does not fit conventions.")
        else: print("Error! Can't find sense sequence in reference.")
    elif sgRNA_strand == 'antisense':
        refSeq_rev = rev_complement(refSeq)
        pos_in_exon = refSeq_rev.find(sgRNA_seq)
        if pos_in_exon >= 0:
            if gene_strand =='plus':
                return refSeq_no[1] - pos_in_exon
            elif gene_strand == 'minus':
                return refSeq_no[0] + pos_in_exon
            else: print("Error! Gene strand info does not fit conventions.")
        else: print("Error! Can't find antisense sequence for", sgRNA_seq ,"in reference.")
    else: print("Error! Strand information not available.")

# Parse seqId so that only exon number remains
def parseName(seqId):
    seqId = int(((seqId.split(' ')[0]).split('_'))[4])
    return seqId

# Parse sequence description to pull out genomic coordinates
def parsePos(seqDesc):
    genPos = ((seqDesc.split(' ')[1]).split(':')[1]).split('-')
    return (int(genPos[0]), int(genPos[1]))

# Find all indexes of a particular base in edit window
def find_Indexes(base, seq_str, indexList):
    last_found = seq_str.rfind(base)
    if last_found == -1:
        return indexList
    indexList.append(last_found)
    new_seqstr = seq_str[:last_found]
    return find_Indexes(base, new_seqstr, indexList)

# Make list of all subsets
def generateSubsets(indexList):
    subsetList = []
    for i in range(1,len(indexList)+1):
        subsetList.extend(combinations(indexList,i))
    return subsetList

# Replace base at specified index to new base (note: original script has an error, fixed here)
def replaceBase(index_tuple, seq_str, new_base):
    for index in index_tuple:
        seq_str = seq_str[:index] + new_base + seq_str[index+1:]
    return seq_str

# Generate all possible combinations of mutations
def generateMutCom(seq_str, subsetList, new_base):
    return [replaceBase(index_tuple, seq_str, new_base) for index_tuple in subsetList]
