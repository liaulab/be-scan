import sys
import re

sys.path.append('../../bigscam/helper') # relative path to helper directory
from genomic import DNA_AA_map, base_editing_key, bases, complements, cas_key
from genomic import rev_complement, complement, protein_to_AAseq, process_PAM

# unit tests for functions in helper/genomic.py

def test_complement(): 
    assert complement(complements, 'ACGT') == 'TGCA'
    assert complement(complements, 'AAAA') == 'TTTT'
    assert complement(complements, '') == ''

def test_rev_complement(): 
    assert rev_complement(complements, 'ACGT') == 'ACGT'
    assert rev_complement(complements, 'AAAA') == 'TTTT'
    assert rev_complement(complements, '') == ''

def test_process_PAM():
    assert process_PAM('NNN') == re.compile('([acgtACGT]{1}[acgtACGT]{1}[acgtACGT]{1})')
    assert process_PAM('AAA') == re.compile('([aA]{1}[aA]{1}[aA]{1})')
    assert process_PAM('RYGT') == re.compile('([aAgG]{1}[cCtT]{1}[gG]{1}[tT]{1})')

def test_protein_to_AAseq(): 
    pass
