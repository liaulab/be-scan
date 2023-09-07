"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: some genomic and transcriptomic variables and functions, 
many of these would be included in bioconda but I'm writing these out to avoid messy dependencies}
"""

import re

# variables

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

base_editing_key = {"CBE": ["C", "T"], 
                    "ABE": ["A", "G"]
                    }

bases = 'ACGT'

complements = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 
                   'a':'t', 't':'a', 'g':'c', 'c':'g'}

cas_key = {'Sp': 'NGG', 'SpG': 'NGN', 'SpRY': 'NNN'}

# functions

def DNA_to_AA(seq): 
    aa_seq = ''
    for i in range(len(seq)//3): 
        codon = seq[(i*3):(i*3)+3] #.replace("T", "U")
        aa_seq += DNA_AA_map[codon]
    return aa_seq

def rev_complement(complements, seq): 
    compl = ''
    for i in range(len(seq)): 
        compl += complements[seq[i]]
    return compl[::-1]

def complement(complements, seq): 
    compl = ''
    for i in range(len(seq)): 
        compl += complements[seq[i]]
    return compl

def protein_to_AAseq(filename): 
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
    PAM = PAM.replace("G", "[gG]{1}")
    PAM = PAM.replace("C", "[cC]{1}")
    PAM = PAM.replace("T", "[tT]{1}")
    PAM = PAM.replace("A", "[aA]{1}")
    PAM = PAM.replace("Y", "[cCtT]{1}")
    PAM = PAM.replace("R", "[aAgG]{1}")
    PAM = PAM.replace("N", "[acgtACGT]{1}")
    return re.compile("({})".format(PAM))
