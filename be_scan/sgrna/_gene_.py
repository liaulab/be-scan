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
# from _genomic_ import complements, rev_complement

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
            introns.append((''.join([base for base in exon[:len(exon)//2] if base.islower()]), 
                            ''.join([base for base in exon[len(exon)//2:] if base.islower()])))

        # CHECK INTRON LENGTHS CONSISTENT ACROSS FASTA FILE #
        assert all([len(x)-len(y)==0 for x, y in introns]), [(len(x), len(y)) for x, y in introns]
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
    # SAVES GUIDE SEQ, INDEX OF FIRST BP IN CHROMOSOME, STARTING FRAME (0,1,2), EXON #
    ### conditional on introns being lowercase and exons being uppercase
    def find_all_guides(self, window, n=23): 
        assert len(self.exons_extra) > 0
        first_exon_utr_len = first_lowercase_length(self.exons_extra[0])
        last_exon_utr_len = first_lowercase_length(self.exons_extra[-1])
        self.n = n
        fwd_guides = []
        prev_frame, prev_ind = 0, 0
        # IDENTIFY ALL N LENGTH SEQUENCES IN ALL EXONS #
        utr_bool = True
        for e, exon_extra in enumerate(self.exons_extra): 
            for i in range(len(exon_extra)-self.n-1): 
                seq = exon_extra[i:i+self.n]
                frame = (i+prev_frame-self.intron_len)%3

                if self.strand == 'plus': ind_chr = self.exons_start[e]+i
                if self.strand == 'minus': ind_chr = self.exons_start[e]-i
                if seq[0].islower(): ind_gene = -1
                else: ind_gene = i-self.intron_len+prev_ind

                # if (e==0 and i < self.intron_len-window[1]-1) or (e==len(self.exons_extra) and i > len(exon_extra)-self.intron_len+1): ind_gene = -2 # OUTSIDE GENE UTR -2
                if seq[window[0]-1].islower(): window_pos1 = -1
                else: window_pos1 = i-self.intron_len+prev_ind+window[0]-1
                if seq[window[1]].islower(): window_pos2 = -1
                else: window_pos2 = i-self.intron_len+prev_ind+window[1]-1

                utr, anti_utr = False, False
                if e == 0 and i < first_exon_utr_len - window[0]+1: 
                    utr = True
                if e == len(self.exons_extra)-1 and i > len(exon_extra) - last_exon_utr_len - window[1]: 
                    utr = True
                if e == 0 and i < first_exon_utr_len + 20 - window[0]+1: 
                    anti_utr = True
                if e == len(self.exons_extra)-1 and i > len(exon_extra) - last_exon_utr_len - 20 + window[1] - window[1]: 
                    anti_utr = True

                calc_pos = i-self.intron_len+prev_ind # USED TO CALC FOR REV GUIDES #
                fwd_guides.append([seq, seq[:self.guide_len], seq[self.guide_len:], 
                                   frame, ind_gene, ind_chr, e, window_pos1, window_pos2, utr, calc_pos, anti_utr])
            prev_frame = (prev_frame+len(exon_extra)-(2*self.intron_len))%3
            prev_ind += len(exon_extra)-(2*self.intron_len)

        # UPDATE INSTANCE VARIABLE #
        self.fwd_guides = [g[1:10] for g in fwd_guides]
        # CALCULATE REV GUIDES FROM FWD GUIDES #
        self.rev_guides = []
        for g in fwd_guides: 
            g_rev = [rev_complement(complements, g[0][3:]), rev_complement(complements, g[0][:3]), (g[3]+1)%3]
            if g_rev[0][0].islower(): g_rev.append(-1)
            else: g_rev.append(g[10]+self.n-1)
            if self.strand=='plus': g_rev.append(g[5]+self.n-1)
            if self.strand=='minus': g_rev.append(g[5]-self.n+3)
            g_rev.append(g[6])

            if    g_rev[0][window[0]-1].islower(): g_rev.append(-1)
            else: g_rev.append((g[10]+self.n-1)-(window[0]-1))
            if    g_rev[0][window[1]].islower(): g_rev.append(-1)
            else: g_rev.append((g[10]+self.n-1)-(window[1]-1))

            g_rev.append(g[11])
            self.rev_guides.append(g_rev)

    # EXTRACT METADATA ABOUT CHROMOSOME POSITION, STRAND #
    def extract_metadata(self): 
        all_lines = self.file_content.split('\n')
        first_line = all_lines[0]
        try: 
            first_line_list = first_line.split()
        except AttributeError: 
            print(f"Error processing {self.filepath} file")

        # EXTRACT CHR NUMBER #
        self.chrID = re.findall(r"(chr(X|Y|\d{1,2}))", first_line_list[1])[0][0]
        # EXTRACT STRAND DIRECTION #
        self.strand = 'plus' if ('strand=+' in first_line_list[4]) else 'minus'
        self.exon_metadata = {}

        header_lines = [line for line in all_lines if (line and line[0] == '>')]
        for i, line in enumerate(header_lines): 
            match = re.search(r"range=(chr[\w\d]+):(\d+)-(\d+)", line)
            
            if match:
                chr_id, start, end = match.groups()
                exon_metadata = {'chromosome':chr_id, 'start':start, 'end':end}
                self.exon_metadata[i] = exon_metadata

def first_lowercase_length(sequence):
    match = re.match(r"[a-z]+", sequence)
    return len(match.group()) if match else 0

def last_lowercase_length_regex(sequence):
    match = re.search(r"[a-z]+$", sequence)
    return len(match.group()) if match else 0
