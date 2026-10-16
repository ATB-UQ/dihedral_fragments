"""Offline grammar tests: canonical form invariances, round trips, and the
pattern-matching assertions that test.py only ran under ``__main__``."""
from itertools import permutations

import pytest

from dihedral_fragments.dihedral_fragment import Dihedral_Fragment
from dihedral_fragments.pattern_matching import re_pattern_matching_for, sql_pattern_matching_for


def _s(left, a, b, right, cycles=None):
    atom_list = (list(left), a, b, list(right)) + ((cycles,) if cycles else ())
    return str(Dihedral_Fragment(atom_list=atom_list))


def test_substituent_order_is_irrelevant():
    ref = _s(['O', 'N', 'C'], 'C', 'C', ['H', 'H', 'C'])
    for left in permutations(['O', 'N', 'C']):
        for right in permutations(['H', 'H', 'C']):
            assert _s(left, 'C', 'C', right) == ref


def test_flipping_the_bond_is_irrelevant():
    assert _s(['O', 'H'], 'C', 'N', ['C', 'H', 'H']) == _s(['C', 'H', 'H'], 'N', 'C', ['O', 'H'])


@pytest.mark.parametrize('text', [
    'C,H,H|C|C|H,H,H',
    'O,N,H|C|C|C,C,H|020,121',
    'C,C,C|C|C|C,C,C|000,101,202',
    'CL,H|C|C|BR,H',
    'C4,H1,H1|C4|C4|H1,H1,H1',
])
def test_canonical_string_is_a_fixed_point(text):
    once = str(Dihedral_Fragment(text))
    assert str(Dihedral_Fragment(once)) == once


def test_ring_cycles_follow_substituent_renumbering():
    # the same ring written with the substituents in two orders
    a = _s(['H', 'H', 'C'], 'C', 'C', ['H', 'C', 'H'], [[2, 0, 1]])
    b = _s(['C', 'H', 'H'], 'C', 'C', ['C', 'H', 'H'], [[0, 0, 0]])
    assert a == b == 'C,H,H|C|C|C,H,H|000'


def test_pattern_matching_asserts_from_test_py_main():
    assert re_pattern_matching_for('Z,%|Z|Z|Z,%')('C,H|C|C|C,H')
    assert not re_pattern_matching_for('Z|Z|Z|Z,%')('C,H|C|C|C,H')
    assert re_pattern_matching_for('Z|Z|Z|Z,%')('C|C|C|C,H')
    assert re_pattern_matching_for('N,J|C|C|J{3}')('N,H|C|C|C,H,H')
    assert re_pattern_matching_for('N,J|C|C|J{2}')('N,H|C|C|C,H')


def test_plain_pattern_is_canonicalised_in_sql():
    assert sql_pattern_matching_for('H,H,C|C|C|H,H,H') == 'dihedral_string="C,H,H|C|C|H,H,H"'


@pytest.mark.xfail(strict=True, reason=(
    'review 2026-10 dihedral_fragments F1: element_valence_for_atom splits "CL1" into '
    '("L", 1), so two-letter elements with a valence digit sort as unknown (999) and the '
    'canonical form depends on input order'))
def test_two_letter_elements_with_valence_are_canonical():
    assert _s(['CL1', 'BR1', 'C4'], 'C4', 'C4', ['H1'] * 3) == _s(['BR1', 'CL1', 'C4'], 'C4', 'C4', ['H1'] * 3)


@pytest.mark.xfail(strict=True, reason=(
    'review 2026-10 dihedral_fragments F2: sql_pattern_matching_for interpolates the user '
    'string into SQL without quoting'))
def test_sql_pattern_does_not_pass_quotes_through():
    assert '"' not in sql_pattern_matching_for('C" OR 1=1 OR "|C|C|H').split('"', 1)[1].rsplit('"', 1)[0]
