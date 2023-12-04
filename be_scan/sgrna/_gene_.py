"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: this class is meant to create a Gene object by preprocessing a .fasta file
This preprocessing includes separating intron and exons, and generating all possible guides of length n}
"""

import re

from be_scan.sgrna._genomic_ import complements
from be_scan.sgrna._genomic_ import rev_complement, complement

# this function is conditional on +-20 bps of introns flanking exons
class GeneForCRISPR(): 
    
    # initialize the .fasta file by reading the file into a string
    def __init__(self, filepath, output_dir=''): 
        self.filepath = filepath
        f = open(filepath, "r")
        self.file_content = f.read()
    
    # function to read in a .fasta file as its exons, separating into just exon and exon +- 20 bps
    # outputs: none
    # effect: saves a list of exons, list of exons +- 20 bps
    # conditional on introns being lowercase and exons being uppercase
    def parse_exons(self): 
        exons_info = []
        exons_extra = []
        i = -1
        # split up file contents into strings and add line by line to exon+-20bps list
        for line in self.file_content.split('\n'): 
            if len(line) > 0 and line[0] == '>': # it title line
                exons_extra.append('')
                exons_info.append(line)
                i += 1
            else: # if bases
                exons_extra[i] += line
        exons = []
        # extract exons with uppercase letters
        for exon in exons_extra: 
            exons.append(''.join([base for base in exon if base.isupper()]))
        # identify the start site of each exon
        exons_start = []
        for info in exons_info: 
            exons_start.append(int(re.findall(r":(\d+)-", info)[0]))
        # change instance variables
        self.exons_extra, self.exons, self.starts = exons_extra, exons, exons_start
    
    # generates all fwd and rev potential guides and their metadata
    # while this is computationally naive, genes are usually not long and can be computationally exhausted easily
    # outputs: none
    # effect: saves guide sequences, the index of the first bp of the guide, 
    #         the frame (0, 1, 2) of the first bp, and the exon # of guide
    def find_all_guides(self, n=23): 
        assert isinstance(n, int)
        self.n = n
        self.fwd_guides = []
        prev_frame = 0
        # identify all n length sequences in exons
        for e, exon_extra in enumerate(self.exons_extra): 
            prev_ind = 0
            for i in range(len(exon_extra)-self.n-1): 
                # the number 20 is due to 20 intron bps assumed to be present
                frame = (i+prev_frame-20)%3
                ind = i+prev_ind+self.starts[e]
                # add to instance variable
                self.fwd_guides.append([exon_extra[i:i+self.n], frame, ind, e])
            prev_frame = (prev_frame+len(exon_extra)-40)%3
        # change instance variables
        self.rev_guides = [[rev_complement(complements, g[0])] + [(g[1]+1)%3] + [g[2]+self.n-1] + [g[3]] for g in self.fwd_guides]
        
    def extract_metadata(self): 
        with open(self.filepath) as f:
            first_line = f.readline()
        try: 
            first_line_list = first_line.split()
        except AttributeError: 
            print(f"Error processing {self.filepath} file")

        # extract chromosome number
        self.chrID = re.findall(r"(chr(X|Y|\d{1,2}))", first_line_list[1])[0][0]
        # extract strand information
        self.strand = 'plus' if ('+' in first_line_list[4]) else 'minus'


