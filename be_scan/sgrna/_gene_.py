"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: this class is meant to create a Gene object by preprocessing a .fasta file
This preprocessing includes separating intron and exons, and generating all possible guides of length n}
"""

import re
import warnings

from pathlib import Path
from be_scan.sgrna._genomic_ import complements, rev_complement

# CLASS OBJECT NEEDED TO PARSE SPECIFIC INPUT INTRON/EXON .FASTA FORMAT #
class GeneForCRISPR(): 
    
    # READING THE .FASTA FILE INTO A STRING #
    def __init__(self, filepath): 
        self.filepath = Path(filepath)
        f = open(self.filepath, "r")
        self.file_content = f.read() # STRING
        self.guide_len = 20

        if len(self.file_content) == 0: 
            warnings.warn('Input file is empty.')
    
    # PARSE SEPARATE EXONS, AND FIND WHERE INTRON/EXON BOUNDARIES ARE #
    # OUTPUT: NONE #
    # SAVES A LIST OF EXONS AND A LIST OF EXONS + INTRON ENDS #
    def parse_exons(self): 
        exons_info, exons_extra = [], []
        i = -1
        # SPLIT FILE CONTENT INTO STRINGS, ADD TO LISTS #
        for line in self.file_content.split('\n'): 
            if len(line) > 0 and line[0] == '>': # TITLE LINE
                exons_extra.append('')
                exons_info.append(line)
                i += 1
            else: # BASES LINE
                exons_extra[i] += line
        # EXTRACT UPPERCASE (EXONS) and LOWERCASE (INTRONS) #
        exons, introns = [], []
        for exon in exons_extra: 
            exons.append(''.join([base for base in exon if base.isupper()]))
            introns.append((''.join([base for base in exon[:len(exon)//2] if base.isupper()]), 
                            ''.join([base for base in exon[len(exon)//2:] if base.isupper()])))

        # CHECK INTRON LENGTHS CONSISTENT ACROSS FASTA FILE #
        assert all([len(x)-len(y)==0 for x, y in introns]), "Make sure all flanking introns are the same length"
        assert all([len(x) == len(y) for x, y in introns]), "Make sure 5' and 3' intron sequences are the same length" ###
        self.intron_len = len(introns[0][0])
        
        # IDENTIFY START POSITION OF EACH EXON #
        exons_start = []
        for info in exons_info: 
            exons_start.append(int(re.findall(r":(\d+)-", info)[0]))
        # SET INSTANCE VARIABLES #
        self.exons_extra, self.exons, self.exons_start = exons_extra, exons, exons_start
    
    # GENERATE ALL FWD AND REV GUIDES AND METADATA #
    # OUTPUT: NONE #
    # effect: saves guide sequences, the index of the first bp of the guide in the gene and in the chromosome, 
    #         the frame (0, 1, 2) of the first bp, and the exon # of guide
    ### conditional on introns being lowercase and exons being uppercase
    def find_all_guides(self, n=23): 
        self.n = n
        fwd_guides = []
        prev_frame = 0
        # IDENTIFY ALL N LENGTH SEQUENCES #
        for e, exon_extra in enumerate(self.exons_extra): 
            for i in range(len(exon_extra)-self.n-1): 
                seq = exon_extra[i:i+self.n]
                frame = (i+prev_frame-self.intron_len)%3
                ind_chr = i+self.exons_start[e]
                fwd_guides.append([seq, seq[:self.guide_len], seq[self.guide_len:], frame, ind_chr, e])
            prev_frame = (prev_frame+len(exon_extra)-(2*self.intron_len))%3

        # UPDATE INSTANCE VARIABLE #
        self.fwd_guides = [g[1:] for g in fwd_guides]
        self.rev_guides = [[rev_complement(complements, g[0][3:]), rev_complement(complements, g[0][:3]), 
                           (g[3]+1)%3, g[5]+self.n-1, g[6]] for g in fwd_guides]

    # EXTRACT METADATA ABOUT CHROMOSOME POSITION, STRAND #
    def extract_metadata(self): 
        first_line = self.file_content.split('\n')[0]
        try: 
            first_line_list = first_line.split()
        except AttributeError: 
            print(f"Error processing {self.filepath} file")

        # EXTRACT CHR NUMBER #
        self.chrID = re.findall(r"(chr(X|Y|\d{1,2}))", first_line_list[1])[0][0]
        # EXTRACT STRAND DIRECTION #
        self.strand = 'plus' if ('strand=+' in first_line_list[4]) else 'minus'
