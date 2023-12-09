
import pytest
from be_scan.sgrna._guideRNA_ import mutation_combos, format_mutation, annotate_mutations, categorize_mutations, calc_editing_window, calc_coding_window, calc_target

def test_make_mutations():
    mutated = mutation_combos("ATCG", edit=["C", "T"], dir="sense")
    assert set(mutated) == {'ATCG', 'ATTG'}

def test_make_mutations_processive():
    mutated = mutation_combos("ATCCG", edit=["C", "T"], dir="sense")
    assert set(mutated) == {'ATCCG', 'ATTCG', 'ATCTG', 'ATTTG'}