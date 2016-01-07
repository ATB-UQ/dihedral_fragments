from copy import deepcopy, copy
from itertools import product
import sys
sys.path.insert(0, '../')
from atb_helpers.iterables import group_by
from itertools import permutations

CHEMICAL_GROUPS = (
    ('carboxylic acid', '%,O|C|O|H'),
    ('ester', '%,O|C|O|C'),
    ('alkyne', '_|C|C|_'),
    ('alkene', '_,_|C|C|_,_'),
    ('alcool I', 'C,H,H|C|O|H'),
    ('alcool II', 'C,C,H|C|O|H'),
    ('alcool III', 'C,C,C|C|O|H'),
    ('amine I', '%|C|N|H,H'),
    ('amine II', '%|C|N|C,H'),
    ('amine III', '%|C|N|C,C'),
    ('monofluoro', 'F,%|C|%|%'),
    ('difluoro', 'F,F,%|C|%|%'),
    ('trifluoro', 'F,F,F|C|%|%'),
    ('chloro', 'CL,%|C|%|%'),
    ('dichloro', 'CL,CL,[^C]*[^L]|C|%|%'),
    ('trichloro', 'CL,CL,CL|C|%|%'),
    ('bromo', 'BR,%|C|%|%'),
    ('dibromo', 'BR,BR,%|C|%|%'),
    ('tribromo', 'BR,BR,BR|C|%|%'),
    ('iodo', 'I,%|C|%|%'),
    ('diiodo', 'I,I,%|C|%|%'),
    ('triiodo', 'I,I,I|C|%|%'),
    ('thiol', '%|C|S|H'),
)

class FragmentDihedral(object):

    def __init__(self, dihedral_string=None, atom_list=None):

        if dihedral_string:
            splitted_string = dihedral_string.split("|")
            self.neighbours_1 = [atom.upper() for atom in splitted_string[0].split(',')]
            self.atom_2 = splitted_string[1].upper()
            self.atom_3 = splitted_string[2].upper()
            self.neighbours_4 = [atom.upper() for atom in splitted_string[3].split(',')]
        else:
            self.neighbours_1, self.atom_2, self.atom_3, self.neighbours_4= atom_list

        canonical_rep = self.__canonical_rep__()

        self.neighbours_1 = canonical_rep.neighbours_1
        self.atom_2 = canonical_rep.atom_2
        self.atom_3 = canonical_rep.atom_3
        self.neighbours_4 = canonical_rep.neighbours_4

    def __str__(self):
        return "{neighbours_1}|{atom_2}|{atom_3}|{neighbours_4}".format(
            neighbours_1=','.join(self.neighbours_1),
            atom_2=self.atom_2,
            atom_3=self.atom_3,
            neighbours_4=','.join(self.neighbours_4),
        )

    def __eq__(self, other):
        return self.__dict__ == other.__dict__
        #return self.__canonical_rep__().__dict__ == other.__canonical_rep__().__dict__

    def __ne__(self, other):
        return not self == other

    def __canonical_rep__(self):
        other = copy(self)
        # Order each neighbour list by alphabetical order
        # WARNING: Keep in mind that maintaining order will be necessary for maintaining sterochemistry information
        other.neighbours_1.sort()
        other.neighbours_4.sort()
        # Compare the two neighbour list and put the first in lexical order on the leftmost side
        if [other.atom_2] + other.neighbours_1 >=  [other.atom_3] + other.neighbours_4:
            other.atom_2, other.atom_3 = other.atom_3, other.atom_2
            other.neighbours_1, other.neighbours_4 = other.neighbours_4, other.neighbours_1
        return other

SQL_SUBSTITUTION_CHARACTERS = ('_', '%')
SQL_FULL_REGEX_CHARACTERS = ('.', '[', ']', '^', '$')
has_substitution_pattern = lambda x: any([y in x for y in SQL_SUBSTITUTION_CHARACTERS])
has_regex_pattern = lambda x: any([y in x for y in SQL_FULL_REGEX_CHARACTERS])

REGEX_START, REGEX_END = ('^', '$')

def escaped_special_regex_characters(patterns):
    return [ (REGEX_START + pattern.replace('%', '[A-Z,]*').replace('|', '\\\\|') + REGEX_END) for pattern in patterns]

def sql_OR(*args):
    return ' '.join(
        ['('] + [' OR '.join(*args)] + [')']
    )

TRUE_OR_FALSE = [False, True]
FALSE = [False]
PUT_SUBSTITUTION_PATTERN_FIRST = True

def sql_pattern_matching_for(pattern):
    assert len([x for x in pattern if has_substitution_pattern(x) ]) <= 5, [x for x in pattern if has_substitution_pattern(pattern) ]

    if not (has_substitution_pattern(pattern) or has_regex_pattern(pattern)):
        return 'dihedral_string="{pattern}"'.format(pattern=pattern)
    else:
        use_full_regex = has_regex_pattern(pattern)

        components = pattern.split('|')
        assert len(components) == 4
        need_to_reverse_inner_atoms = (components[1] != components[2])

        left_permutations = permutations(range(len(group_by(components[0].split(','), key=lambda x:x))))
        right_permutation = permutations(range(len(group_by(components[3].split(','), key=lambda x:x))))

        patterns = [
            correct_pattern(pattern, should_reverse=should_reverse, left_permutation=left_permutation, right_permutation=right_permutation)
            for (should_reverse, left_permutation, right_permutation) in product((TRUE_OR_FALSE if need_to_reverse_inner_atoms else FALSE), left_permutations, right_permutation)
        ]


        if use_full_regex:
            patterns = escaped_special_regex_characters(patterns)

        return sql_OR(
            [
                'dihedral_string {sql_operator} "{match_pattern}"'.format(
                    sql_operator=('LIKE' if not use_full_regex else 'REGEXP'),
                    match_pattern=match_pattern,
                )
                for match_pattern in patterns
            ]
        )

def correct_pattern(pattern, should_reverse=False, left_permutation=(), right_permutation=()):
    return '|'.join(
        (reversed if should_reverse else lambda x:x)(
            [sorted_components(i, x, left_permutation=left_permutation, right_permutation=right_permutation) for (i, x) in enumerate(pattern.split('|'))]
        )
    )

def sorted_components(component_index, component, left_permutation=(), right_permutation=()):
    if component_index in (0,3):
        return ','.join(
            sorted_components_list(
                component.split(','),
                permutation=left_permutation if component_index == 0 else right_permutation,
            )
        )
    else:
        return component

def sorted_components_list(component_list, permutation=()):

    sorting_dict = dict(
        zip(
            sorted(group_by(component_list, lambda x:x).keys()),
            permutation,
        )
    )

    return sorted(
        component_list,
        key=lambda x: sorting_dict[x],
    )

if __name__ == "__main__" :
    for pattern in ('C|N||C|C', 'C,A,B,D|N|N|H,_,C', 'CL,CL,[^C]*[^L]*|C|%|%', 'CL,CL,[^C]*[^L]*|C|C|%'):
        print pattern
        print sql_pattern_matching_for(pattern)
        print
    exit()

    dihedral_1 = FragmentDihedral("C,C,H|C|C|C,H,H")
    print dihedral_1
    dihedral_2 = FragmentDihedral("C,H,H|C|C|C,C,H")
    print dihedral_2
    print dihedral_1 == dihedral_2
    dihedral_3 = FragmentDihedral(atom_list=(['H', 'H'], 'C', 'C', ['Cl', 'Cl']))
    print dihedral_3
    print dihedral_3 == dihedral_2


