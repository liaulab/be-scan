"""
Author: Calvin XiaoYang Hu
Date: 230907

{Description: suit of unit tests for testing functions in be_scan/be_scan/helper/genomic}
"""

import re
import pytest

from be_scan.sgrna._genomic_ import bases, complements, cas_key
from be_scan.sgrna._genomic_ import rev_complement, complement, protein_to_AAseq, process_PAM, DNA_to_AA

# unit tests for functions in helper/genomic.py

def test_complement_pos(): # positive control test
    assert complement(complements, 'ACGT') == 'TGCA'
    assert complement(complements, 'AAAA') == 'TTTT'
    assert complement(complements, '') == ''

def test_complement_neg(): # negative control test
    with pytest.raises(AssertionError): # sequence has non ACGT letters
        assert complement(complements, 'ABCD') == 'TBGD'
    with pytest.raises(AssertionError): # sequence is not a string
        assert complement(complements, 1) == ''

def test_rev_complement_pos(): # positive control test
    assert rev_complement(complements, 'ACGT') == 'ACGT'
    assert rev_complement(complements, 'AAAA') == 'TTTT'
    assert rev_complement(complements, '') == ''

def test_rev_complement_neg(): # negative control test
    with pytest.raises(AssertionError): # sequence has non ACGT letters
        assert rev_complement(complements, 'ABCD') == 'DGBT'
    with pytest.raises(AssertionError): # sequence is not a string
        assert rev_complement(complements, 1) == ''

def test_process_PAM_pos(): # positive control test
    assert process_PAM('NNN') == re.compile('([acgtACGT]{1}[acgtACGT]{1}[acgtACGT]{1})')
    assert process_PAM('AAA') == re.compile('([aA]{1}[aA]{1}[aA]{1})')
    assert process_PAM('RYGT') == re.compile('([aAgG]{1}[cCtT]{1}[gG]{1}[tT]{1})')
    assert process_PAM('ryGt') == re.compile('([aAgG]{1}[cCtT]{1}[gG]{1}[tT]{1})')

def test_process_PAM_neg(): # negative control test
    with pytest.raises(AssertionError): # sequence has non ACGT letters
        assert process_PAM('abcd') == re.compile('([acgtACGT]{1}[acgtACGT]{1}[acgtACGT]{1})')
    with pytest.raises(AssertionError): # sequence is not a string
        assert process_PAM(1) == re.compile('([acgtACGT]{1}[acgtACGT]{1}[acgtACGT]{1})')

def test_protein_to_AAseq_pos(): # positive control test
    protein_seq = protein_to_AAseq(filename='tests/test_data/sgrna_data/P10275.fasta')
    assert len(protein_seq) == 921

def test_protein_to_AAseq_neg(): # negative control test
    protein_seq = protein_to_AAseq(filename='tests/test_data/sgrna_data/P10275_changed.fasta')
    assert len(protein_seq) == 360

def test_DNA_to_AA_pos(): # positive control test
    assert DNA_to_AA('TTTAAACCCGGG') == 'FKPG'
    assert DNA_to_AA('') == ''

def test_DNA_to_AA_neg(): # negative control test
    with pytest.raises(AssertionError): # sequence has non ACGT letters
        assert DNA_to_AA('TTTAAACCCabc') == 'FKPG'
    with pytest.raises(AssertionError): # sequence is not a string
        assert DNA_to_AA(1) == 'FKPG'
    with pytest.raises(AssertionError): # sequence is not right length
        assert DNA_to_AA('cgttt') == 'FKPG'
