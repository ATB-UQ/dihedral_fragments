from copy import deepcopy, copy
from itertools import product
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

def sql_OR(*args):
    return ' '.join(
        ['('] + [' OR '.join(*args)] + [')']
    )

TRUE_OF_FALSE = [False, True]
PUT_SUBSTITUTION_PATTERN_FIRST = True

def sql_pattern_matching_for(pattern):
    assert len([x for x in pattern if has_substitution_pattern(x) ]) <= 3, [x for x in pattern if has_substitution_pattern(pattern) ]

    if not has_substitution_pattern(pattern):
        return 'dihedral_string="{pattern}"'.format(pattern=pattern)
    else:
        use_full_regex = has_regex_pattern(pattern)
        patterns = [
            correct_pattern(pattern, should_reverse=should_reverse, put_substitution_pattern_first=put_substitution_pattern_first)
            for (should_reverse, put_substitution_pattern_first) in product(TRUE_OF_FALSE, TRUE_OF_FALSE)
        ]
        return sql_OR(
            [
                'dihedral_string {sql_operator} "{match_pattern}"'.format(
                    sql_operator=('LIKE' if not use_full_regex else 'REGEXP'),
                    match_pattern=match_pattern,
                )
                for match_pattern in patterns
            ]
        )

def correct_pattern(pattern, should_reverse=False, put_substitution_pattern_first=PUT_SUBSTITUTION_PATTERN_FIRST):
    return '|'.join(
        (reversed if should_reverse else lambda x:x)(
            [sorted_components(i, x, put_substitution_pattern_first=put_substitution_pattern_first) for (i, x) in enumerate(pattern.split('|'))]
        )
    )

def sorted_components(component_index, component, put_substitution_pattern_first=PUT_SUBSTITUTION_PATTERN_FIRST):
    if component_index in (0,3):
        return ','.join(
            sorted_components_list(
                component.split(','),
                put_substitution_pattern_first=put_substitution_pattern_first,
            )
        )
    else:
        return component

def sorted_components_list(component_list, put_substitution_pattern_first=PUT_SUBSTITUTION_PATTERN_FIRST):
    sort_first = lambda x: (has_regex_pattern(x) or has_substitution_pattern(x))

    if put_substitution_pattern_first:
        final_sort_first = lambda x: not sort_first(x)
    else:
        final_sort_first = lambda x: sort_first(x)

    return sorted(
        component_list,
        key=lambda x: (final_sort_first(x), x)
    )

if __name__ == "__main__" :
    print sql_pattern_matching_for('C|N||C|C')
    print sql_pattern_matching_for('C,A|N|N|H,_,C')
    exit()

    dihedral_1 = FragmentDihedral("C,C,H|C|C|C,H,H")
    print dihedral_1
    dihedral_2 = FragmentDihedral("C,H,H|C|C|C,C,H")
    print dihedral_2
    print dihedral_1 == dihedral_2
    dihedral_3 = FragmentDihedral(atom_list=(['H', 'H'], 'C', 'C', ['Cl', 'Cl']))
    print dihedral_3
    print dihedral_3 == dihedral_2

