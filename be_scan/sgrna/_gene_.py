"""
Author: Calvin XiaoYang Hu
Date: 230906

{Description: this class is meant to create a Gene object by preprocessing a .fasta file
This preprocessing includes separating intron and exons, and generating all possible guides of length n}
"""

import re

from pathlib import Path
from be_scan.sgrna._genomic_ import complements, rev_complement

# this function is conditional on +-20 bps of introns flanking exons
class GeneForCRISPR(): 
    
    # initialize the .fasta file by reading the file into a string
    def __init__(self, filepath, output_dir=''): 
        self.filepath = Path(filepath)
        f = open(self.filepath, "r")
        self.file_content = f.read()
    
    # function to read in a .fasta file as its exons, separating into just exon and exon +- 20 bps
    # outputs: none
    # effect: saves a list of exons, list of exons +- 20 bps
    # conditional on introns being lowercase and exons being uppercase
    def parse_exons(self): 
        exons_info = []
        exons_extra = []
        intron_lens = []
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
            intron_lens.append(sum(1 for c in exon if c.islower()))

        # check intron lengths are consistent across fasta file
        # lengths of introns are consistent across exons
        assert len(set(intron_lens)) == 1, "Make sure all flanking introns are the same length"
        # lengths of introns 5' and 3' are consistent
        assert intron_lens[0] % 2 == 0, "Make sure 5' and 3' intron sequences are the same length"
        self.intron_len = int(intron_lens[0]/2)
        # identify the start site of each exon
        exons_start = []
        for info in exons_info: 
            exons_start.append(int(re.findall(r":(\d+)-", info)[0]))
        # change instance variables
        self.exons_extra, self.exons, self.starts = exons_extra, exons, exons_start
    
    # generates all fwd and rev potential guides and their metadata
    # while this is computationally naive, genes are usually not long and can be computationally exhausted easily
    # outputs: none
    # effect: saves guide sequences, the index of the first bp of the guide in the gene and in the chromosome, 
    #         the frame (0, 1, 2) of the first bp, and the exon # of guide
    def find_all_guides(self, n=23): 
        assert isinstance(n, int)
        self.n = n
        fwd_guides = []
        prev_frame, prev_ind = 0, 0
        # identify all n length sequences in exons
        for e, exon_extra in enumerate(self.exons_extra): 
            for i in range(len(exon_extra)-self.n-1): 
                # the number 20 is due to 20 intron bps assumed to be present
                frame = (i+prev_frame-self.intron_len)%3
                ind_gene = i-self.intron_len+prev_ind
                ind_chr = i+self.starts[e]
                # add to instance variable
                seq = exon_extra[i:i+self.n]
                fwd_guides.append([seq, seq[:20], seq[20:], frame, ind_gene, ind_chr, e])
            prev_frame = (prev_frame+len(exon_extra)-(2*self.intron_len))%3
            prev_ind += len(exon_extra)-(2*self.intron_len)
        # change instance variables
        self.fwd_guides = [g[1:] for g in fwd_guides]
        self.rev_guides = [[rev_complement(complements, g[0][3:]), rev_complement(complements, g[0][:3]), (g[3]+1)%3, g[4]+self.n-1, g[5]+self.n-1, g[6]] for g in fwd_guides]

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
