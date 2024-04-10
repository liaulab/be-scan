"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: some genomic and transcriptomic variables and functions, 
many of these would be included in bioconda etc but I'm writing these out to avoid messy dependencies}
"""

import re
from pathlib import Path

# VARIABLES

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

bases = 'ACGT'
complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 
                   'a':'t', 't':'a', 'g':'c', 'c':'g'}



### add functionality for dual base editos and for R or Y

# name of cas protein with their associated PAM
cas_key = {'Sp': 'NGG', 
           'SpG': 'NGN', 
           'SpRY': 'NNN', 
           'SpRY_highefficiency': 'NRN', 
           'SpRY_lowefficiency': 'NYN'
           }

# FUNCTIONS

# translate DNA sequence to amino acid sequence
def DNA_to_AA(seq, upper=True): 
    assert isinstance(seq, str)
    assert all(c in 'acgtACGT' for c in seq)
    assert len(seq) % 3 == 0
    aa_seq = ''
    for i in range(len(seq)//3): 
        codon = seq[(i*3):(i*3)+3] #.replace("T", "U")
        # if upper is False, we do not translate lowercase letters
        if upper: 
            codon = codon.upper()
        if not codon.isupper(): 
            aa_seq += '_'
            continue
        aa_seq += DNA_AA_map[codon]
    return aa_seq

# find the reverse complement of a DNA sequence
def rev_complement(complements, seq): 
    assert isinstance(seq, str)
    assert all(c in 'acgtACGT' for c in seq)
    compl = ''
    for i in range(len(seq)): 
        compl += complements[seq[i]]
    return compl[::-1]

# find the complement of a DNA sequence
def complement(complements, seq): 
    assert isinstance(seq, str)
    assert all(c in 'acgtACGT' for c in seq)
    compl = ''
    for i in range(len(seq)): 
        compl += complements[seq[i]]
    return compl

# take in a protein .fasta file and extract the protein sequence
def protein_to_AAseq(filename): 
    filepath = Path(filename)
    f = open(filename, "r")
    file_content = f.read().split('\n')
    seq = ''.join(file_content[1:])
    
    dic = {}
    for i in range(len(seq)):
        dic[i+1] = seq[i]
    dic[len(seq)+1] = '.'
    return dic

# function to change a PAM sequence into a regex sequence
def process_PAM(PAM): 
    assert isinstance(PAM, str)
    assert all(c in 'acgtrynACGTRYN' for c in PAM)
    ### make sure only nucleic acid letters are allowed, and to uppercase
    PAM = PAM.upper()
    PAM = PAM.replace("G", "[gG]{1}")
    PAM = PAM.replace("C", "[cC]{1}")
    PAM = PAM.replace("T", "[tT]{1}")
    PAM = PAM.replace("A", "[aA]{1}")
    PAM = PAM.replace("Y", "[cCtT]{1}")
    PAM = PAM.replace("R", "[aAgG]{1}")
    PAM = PAM.replace("N", "[acgtACGT]{1}")
    return re.compile("({})".format(PAM))
