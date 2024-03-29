# SUMMARY: 

# the idea of this workflow is to first create a 
# gene object that has all possible guides listed (_gene_.py)

# then based on CRISPR criteria filter out guides that won't work 
#    (_BE_guides_.py identify_guides)

# then annotate the remaining guides and output them
#    (_BE_guides_.py annotate_guides)

from be_scan.sgrna._genomic_ import DNA_to_AA, rev_complement, complement, protein_to_AAseq, process_PAM, make_mutations
from be_scan.sgrna._BE_guides_ import identify_BE_guides, annotate_BE_guides
# from be_scan.sgrna.findall_be import 
