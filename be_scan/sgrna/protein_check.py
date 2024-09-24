from Bio import SeqIO
from Bio.Seq import Seq
import os

def translate(gene_filepath):
    """[Summary]
    Generate protein sequence from gene sequence

    Parameters
    ------------
    gene_filepath: path
        The file with the gene sequence

    Returns
    ------------
    Seq
        The protein sequence
    """

    exonSeq = ""
    with open(gene_filepath, "r") as gene_file:
        #For each exon in the gene file, add the exon to exonSeq and translate full
        for record in SeqIO.parse(gene_file, "fasta"):
            exon = ''.join([c for c in record if c.isupper()])
            exonSeq += exon

    exonSeq = Seq(exonSeq).translate()

    return exonSeq

# Helper function for checking if given gene and protein align
def protein_check(gene_filepath, protein_filepath):
    """[Summary]
    Check if the gene and protein sequences align

    Parameters
    ------------
    gene_filepath: path
        The file with the gene sequence
    protein_filepath: path
        The file with the protein .fasta sequence

    Returns
    ------------##
    bool
        True if the gene and protein sequences align
    """

    #check if files exist
    if(not os.path.exists(gene_filepath)):
        raise FileNotFoundError("Gene file not found")

    #Translate gene sequence
    prot_translate = translate(gene_filepath)
    prot_seq = ""

    with open(protein_filepath, "r") as protein_file:
        for record in SeqIO.parse(protein_file, "fasta"):
            prot_seq = record.seq

    #hardcode stripped final stop codon :C look for a better way to do this
    if prot_translate[-1] == "*":
        prot_translate = prot_translate[:-1]

    return prot_translate == prot_seq


#TESTING
if __name__ == '__main__':
    #Changing filepath based on file structure of current repo
    gene_filepath='../../tests/test_data/sgrna/230408_AR_Input.fasta'
    protein_filepath='../../tests/test_data/sgrna/P10275.fasta'


    #print(translate(gene_filepath))
    print(protein_check(gene_filepath, protein_filepath))
