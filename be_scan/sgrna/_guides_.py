"""
Author: Calvin XiaoYang Hu
Date: 231102

{Description: helper functions for processing guides and library data}
"""

import re

# evaluates if guide has PAM and has a residue in window
# returns TRUE or FALSE
def filter_guide(g, PAM_regex, PAM, edit, window): 
    return PAM_regex.match(g[0][-len(PAM):]) and edit[0] in g[0][window[0]-1:window[1]]

# delete duplicate guides since duplicate guides are difficult to deconvolute
# input: a list of guides
# output a shorter list of guides with no duplicates
def filter_repeats(results): 
    seqs = [g[0] for g in results]
    dup_seqs = set([x for x in seqs if seqs.count(x) > 1])
    results = [g.copy() for g in results if g[0] not in dup_seqs]
    return results
