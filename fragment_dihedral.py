from copy import deepcopy, copy
from itertools import product
from itertools import permutations
from jinja2 import Template
import re
from collections import namedtuple

from atb_helpers.iterables import group_by
from pipeline.pipelineHelperFunctions import ELEMENT_NUMBERS

DEBUG = False

def print_if_DEBUG(something):
    if DEBUG:
        print something

NEIGHBOUR_SEPARATOR = ','

def join_neighbours(neighbours):
    return NEIGHBOUR_SEPARATOR.join(neighbours)

def split_neighbour_str(neighbour_str):
    return neighbour_str.split(NEIGHBOUR_SEPARATOR)

GROUP_SEPARATOR = '|'

def join_groups(groups):
    return GROUP_SEPARATOR.join(groups)

def split_group_str(group_str):
    return group_str.split(GROUP_SEPARATOR)

class Dihedral_Fragment_Matching_Pattern(object):
    def __init__(self, pattern_string=None):
        splitted_string = split_group_str(pattern)
        self.neighbours_1 = [atom.upper() for atom in split_neighbour_str(splitted_string[0])]
        self.atom_2 = splitted_string[1].upper()
        self.atom_3 = splitted_string[2].upper()
        self.neighbours_4 = [atom.upper() for atom in split_neighbour_str(splitted_string[3])]

def has_regex_pattern(pattern):
    assert ('$' not in pattern) and ('^' not in pattern)
    return any([y in pattern for y in SQL_FULL_REGEX_CHARACTERS])

REGEX_START_ANCHOR, REGEX_END_ANCHOR, REGEX_OR_OPERATOR = ('^', '$', '|')

BACKSLASH = '\\'
COMMA = ','
ESCAPED_COMMA = ';'
UNESCAPE_COMMA = lambda x: x.replace(ESCAPED_COMMA, COMMA)
ATOM_CHARACTERS, MAX_ATOM_CHARACTERS = 'A-Z', 2
VALENCE_CHARACTERS, MAX_VALENCE_CHARACTERS = '0-9', 1
ONE_ATOM = '[' + ATOM_CHARACTERS + VALENCE_CHARACTERS + ']' + '{1,' + str(MAX_ATOM_CHARACTERS + MAX_VALENCE_CHARACTERS) + '}'
ONE_NUMBER = '[0-9]'
ANY_NUMBER_OF_ATOMS = '[' + ATOM_CHARACTERS + VALENCE_CHARACTERS + ESCAPED_COMMA + ']*'

def REGEX_NOT(pattern):
    return '[^(' + str(pattern) + ')]'

def REGEX_GROUP(index):
    return BACKSLASH + str(index)

def REGEX_OR(*args):
    return '(' + REGEX_OR_OPERATOR.join(args) + ')'

def REGEX_ESCAPE(pattern, flavour='sql'):
    if flavour == 'sql':
        return BACKSLASH + str(pattern)
        return BACKSLASH + BACKSLASH + str(pattern)
    elif flavour == 're':
        return '[' + str(pattern) + ']'
    else:
        raise Exception()

def REGEX_AT_LEAST(pattern, escape_plus=True):
    return str(pattern) + (BACKSLASH if escape_plus else '') + '+'

def FORMAT_ESCAPED(pattern):
    return pattern.replace('{', '{{').replace('}', '}}')

def FORMAT_UNESCAPED(pattern):
    return pattern.replace('{{', '{').replace('}}', '}')

def NOT(pattern):
    return '!' + str(pattern)

def CAPTURE(pattern):
    return '(' + str(pattern) + ')'

def ESCAPE(pattern):
    return BACKSLASH + str(pattern)

def on_desc_number_electron_then_desc_valence(atom):
    upper_atom = atom.upper()
    match = re.search(
        CAPTURE('[' + ATOM_CHARACTERS + ']') + CAPTURE('[' + VALENCE_CHARACTERS + ']'),
        upper_atom,
    )
    if match:
        element, valence = match.groups()
        valence = int(valence)
    else:
        element, valence = upper_atom, 1
    ASC = lambda x: x
    try:
        return (
            ASC(ELEMENT_NUMBERS[element]),
            ASC(valence),
        )
    except KeyError:
        return (
            999,
        )

class FragmentDihedral(object):
    HIGHEST_ATOMS_FIRST = dict(
        key=on_desc_number_electron_then_desc_valence,
        reverse=True,
    )

    def __init__(self, dihedral_string=None, atom_list=None):
        assert dihedral_string or atom_list, dict(dihedral_string=dihedral_string, atom_list=atom_list)

        if dihedral_string:
            splitted_string = split_group_str(dihedral_string)
            self.neighbours_1 = [atom.upper() for atom in split_neighbour_str(splitted_string[0])]
            self.atom_2 = splitted_string[1].upper()
            self.atom_3 = splitted_string[2].upper()
            self.neighbours_4 = [atom.upper() for atom in split_neighbour_str(splitted_string[3])]
        else:
            self.neighbours_1, self.atom_2, self.atom_3, self.neighbours_4 = atom_list

        canonical_rep = self.__canonical_rep__()

        self.neighbours_1 = canonical_rep.neighbours_1
        self.atom_2 = canonical_rep.atom_2
        self.atom_3 = canonical_rep.atom_3
        self.neighbours_4 = canonical_rep.neighbours_4

    def __str__(self):
        return "{neighbours_1}|{atom_2}|{atom_3}|{neighbours_4}".format(
            neighbours_1=join_neighbours(self.neighbours_1),
            atom_2=self.atom_2,
            atom_3=self.atom_3,
            neighbours_4=join_neighbours(self.neighbours_4),
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
        other.neighbours_1.sort(**FragmentDihedral.HIGHEST_ATOMS_FIRST)
        other.neighbours_4.sort(**FragmentDihedral.HIGHEST_ATOMS_FIRST)

        # Compare the two central atoms and put the heavier one on the left
        if on_desc_number_electron_then_desc_valence(other.atom_2) <  on_desc_number_electron_then_desc_valence(other.atom_3):
            should_reverse = True
        elif on_desc_number_electron_then_desc_valence(other.atom_2) == on_desc_number_electron_then_desc_valence(other.atom_3):
            should_reverse = False

            if len(self.neighbours_1) < len(self.neighbours_4):
                should_reverse = True
            elif len(self.neighbours_1) == len(self.neighbours_4):
                # If identical central atoms, and same number of neighbours on both ends, try to resolve ambiguity one neighbour at a time
                for (neighbour_1, neighbour_4) in zip(self.neighbours_1, self.neighbours_4):
                    print neighbour_1
                    print on_desc_number_electron_then_desc_valence(neighbour_1)
                    print neighbour_4
                    print on_desc_number_electron_then_desc_valence(neighbour_4)
                    print on_desc_number_electron_then_desc_valence(neighbour_1) < on_desc_number_electron_then_desc_valence(neighbour_4)
                    if on_desc_number_electron_then_desc_valence(neighbour_1) < on_desc_number_electron_then_desc_valence(neighbour_4):
                        should_reverse = True
                        break
                    elif on_desc_number_electron_then_desc_valence(neighbour_1) == on_desc_number_electron_then_desc_valence(neighbour_4):
                        pass
                    else:
                        break
            else:
                should_reverse = False
        else:
            should_reverse = False

        if should_reverse:
            other.reverse_dihedral()
        return other

    def reverse_dihedral(self):
        self.atom_2, self.atom_3 = self.atom_3, self.atom_2
        self.neighbours_1, self.neighbours_4 = self.neighbours_4, self.neighbours_1

Operator_Pattern = namedtuple('Operator_Pattern', 'pattern, replacement, substitution_type')

def exactly_N_times_operator(what=ONE_ATOM, how_many_times=ONE_NUMBER):
    return CAPTURE(what) + ESCAPE('{') + CAPTURE(how_many_times) + ESCAPE('}')

def N_to_M_times_operator(what=ONE_ATOM, N=ONE_NUMBER, M=ONE_NUMBER):
    return CAPTURE(ONE_ATOM) + ESCAPE('{') + CAPTURE(ONE_NUMBER) + ESCAPE('-') + CAPTURE(ONE_NUMBER) + ESCAPE('}')

OPERATORS = (
    (   # !A
        NOT(CAPTURE(ONE_ATOM)),
        REGEX_NOT( REGEX_GROUP(1) ),
        're',
    ),
    (   # A+
        REGEX_AT_LEAST(CAPTURE(ONE_ATOM), escape_plus=True),
        REGEX_AT_LEAST(REGEX_OR( REGEX_GROUP(1), ESCAPED_COMMA), escape_plus=False),
        're',
    ),
    (   # A{X}
        exactly_N_times_operator,
        lambda x, z: FORMAT_ESCAPED( REGEX_OR(x, ESCAPED_COMMA) + '{' + str(2 * int(z) - 1) + '}' ),
        're_map_groups',
    ),
    (   # A{X,Y}
        N_to_M_times_operator,
        lambda x, y, z: FORMAT_ESCAPED( REGEX_OR(x, ESCAPED_COMMA) + '{' + str(int(y) + 1) + ESCAPED_COMMA + str(2 * int(z) - 1) + '}' ),
        're_map_groups',
    ),
)

OPERATORS = map(
    lambda x: Operator_Pattern(*x),
    OPERATORS,
)

def regex_filters(flavour='sql'):
    return [
        # Order matters here !
        Operator_Pattern(
            '|',
            REGEX_ESCAPE('|', flavour=flavour),
            'str',
        ),
    ] + OPERATORS + [
        Operator_Pattern('%', ANY_NUMBER_OF_ATOMS, 'str'),
    ]

OPERATOR_STRIPPER = (
    lambda x: re.find('')
)

ATOM_CATEGORIES = {
    'J': REGEX_OR('C', 'H'),
    'X': REGEX_OR('F', 'I', 'BR', 'CL'),
    'Y': REGEX_OR('N', 'O', 'S'),
    'Z': FORMAT_ESCAPED(ONE_ATOM),
}

SQL_SUBSTITUTION_CHARACTERS = ('_', '%')
SQL_FULL_REGEX_CHARACTERS = ('{', '}', '!', '+', '-') + tuple(ATOM_CATEGORIES.keys())
has_substitution_pattern = lambda x: any([y in x for y in SQL_SUBSTITUTION_CHARACTERS])

def substitute_atom_pattern(atom_pattern):
    def cleaned_atom_pattern(x):
        return re.sub('[\[\]\(\)\^\|\{\}\+\-0-9,]', '', x) # This line prevents atoms from ever having numbers in their name

    for cleaned_atom_pattern in (atom_pattern, cleaned_atom_pattern(atom_pattern)):
        if cleaned_atom_pattern in ATOM_CATEGORIES:
            return re.sub(cleaned_atom_pattern, ATOM_CATEGORIES[cleaned_atom_pattern], atom_pattern)
    return atom_pattern

def substitute_atoms_in_pattern(pattern, flavour='sql'):
    if flavour == 'sql':
        BAR = REGEX_ESCAPE('|')
    elif flavour == 're':
        BAR = '|'
    return BAR.join(
        [
            join_neighbours([ substitute_atom_pattern(atom) for atom in split_on_atoms(group)])
            for group in pattern.split(GROUP_SEPARATOR)
        ]
    )

def split_on_atoms(group):
    atom_patterns = split_neighbour_str(group)

    return atom_patterns

def apply_regex_filters(string, debug=False, flavour='sql'):
    for (pattern, replacement, substitution_type) in regex_filters(flavour):
        if debug:
            old_string = string

        if substitution_type == 'str':
            string = string.replace(pattern, replacement)
        elif substitution_type == 're':
            string = re.sub(pattern, replacement, string)
        elif substitution_type == 're_map_groups':
            general_pattern = pattern()

            matches = re.findall(general_pattern, string)
            if matches:
                if debug:
                    print matches
                for match in matches:
                    tailored_pattern = pattern(*match)
                    mapped_replacement = replacement(*match)
                    if debug:
                        print tailored_pattern, match
                        print mapped_replacement
                    string = re.sub(tailored_pattern, mapped_replacement, string)
        else:
            raise Exception('Wrong substitution_type')

        if debug:
            print [pattern, replacement, old_string], string
    print_if_DEBUG(string)
    return UNESCAPE_COMMA(substitute_atoms_in_pattern(string, flavour=flavour))

def escaped_special_regex_characters(patterns, flavour='sql'):
    return [ (REGEX_START_ANCHOR + apply_regex_filters(pattern, flavour=flavour)  + REGEX_END_ANCHOR) for pattern in patterns]

def sql_OR(*args):
    return ' '.join(
        ['('] + [' OR '.join(*args)] + [')']
    )

TRUE_OR_FALSE = [False, True]
FALSE = [False]
PUT_SUBSTITUTION_PATTERN_FIRST = True

LEFT_GROUP_INDEX, LEFT_ATOM_INDEX, RIGHT_ATOM_INDEX, RIGHT_GROUP_INDEX = (0, 1, 2, 3)

def re_patterns(pattern, full_regex=False, flavour='sql', debug=False, metadata=None):
    components = split_group_str(pattern)
    assert len(components) == 4
    need_to_reverse_inner_atoms = (components[LEFT_ATOM_INDEX] == components[RIGHT_ATOM_INDEX]) or any([x in ATOM_CATEGORIES.keys() for x in (components[LEFT_ATOM_INDEX], components[RIGHT_ATOM_INDEX])])

    if debug:
        print
        print metadata
        print 'need_to_reverse_inner_atoms: {0}'.format(need_to_reverse_inner_atoms)
        print components[LEFT_ATOM_INDEX]
        print components[RIGHT_ATOM_INDEX]

    left_neighbour_groups, right_neighbour_groups = map(
        lambda index: split_neighbour_str(components[index]),
        (LEFT_GROUP_INDEX, RIGHT_GROUP_INDEX),
    )

    left_permutations, right_permutations = map(
        lambda group: permutations(range(len(group_by(group, key=lambda x:x)))),
        (left_neighbour_groups, right_neighbour_groups),
    )

    patterns = [
        correct_pattern(
            pattern,
            should_reverse=should_reverse,
            left_permutation=left_permutation,
            right_permutation=right_permutation,
        )
        for (should_reverse, left_permutation, right_permutation) in product(
            (TRUE_OR_FALSE if need_to_reverse_inner_atoms else FALSE),
            left_permutations,
            right_permutations,
        )
    ]

    if full_regex:
        patterns = escaped_special_regex_characters(patterns, flavour=flavour)

    return patterns

def re_pattern_matching_for(pattern, debug=False, metadata=None):
    if not (has_substitution_pattern(pattern) or has_regex_pattern(pattern)):
        return lambda test_string: test_string == str(FragmentDihedral(pattern))
    else:
        patterns = [FORMAT_UNESCAPED(re_pattern) for re_pattern in re_patterns(pattern, full_regex=True, flavour='re', debug=debug, metadata=metadata)]

        def match_pattern_to(test_string):
            if debug:
                print metadata
                print pattern
                print patterns
                print test_string
                print any([re.search(match_pattern, test_string) for match_pattern in patterns])
                print
            return any([re.search(match_pattern, test_string) for match_pattern in patterns])

        return match_pattern_to

def sql_pattern_matching_for(pattern):
    if not (has_substitution_pattern(pattern) or has_regex_pattern(pattern)):
        return 'dihedral_string="{pattern}"'.format(pattern=FragmentDihedral(pattern))
    else:
        use_full_regex = has_regex_pattern(pattern)
        patterns = re_patterns(pattern, full_regex=use_full_regex, flavour='sql')
        #patterns = [FORMAT_UNESCAPED(re_pattern) for re_pattern in re_patterns(pattern, full_regex=True, flavour='sql')]

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
    return join_groups(
        (reversed if should_reverse else lambda x:x)(
            [sorted_components(i, x, left_permutation=left_permutation, right_permutation=right_permutation) for (i, x) in enumerate(split_group_str(pattern))]
        )
    )

def sorted_components(component_index, component, left_permutation=(), right_permutation=()):
    if component_index in (0,3):
        return join_neighbours(
            sorted_components_list(
                split_neighbour_str(component),
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

SYNTAX_HELP = Template('''
<h5>Syntax Help</h5>

<p class='help block'>
  <span style='font-weight:bold'>Atom categories</span>
  <ul>
    {% for (code, matches) in ATOM_CATEGORIES %}
      <li><code>{{ code }}</code>Any (single) atom in {{ matches }}</li>
    {% endfor %}
  </ul>

  <span style='font-weight:bold'>Operators</span>
  <ul>
    <li><code>!A</code>One single atom not of type <code>A</code>. Ex: <code>!CL|C|C|%</code></li>
    <li><code>A+</code>One or more atoms of type <code>A</code>. Ex: <code>CL+|C|C|%</code></li>
    <li><code>!A+</code>One or more atoms not of type <code>A</code>. Ex: <code>!CL+|C|C|%</code></li>
    <li><code>A{X}</code> Exactly <code>X</code> atoms of type <code>A</code>. Ex: <code>CL{3}|C|C|%</code></li>
    <li><code>A{X-Y}</code> From <code>X</code> to <code>Y</code> atoms of type <code>A</code>. Ex: <code>CL{2-3}|C|C|%</code></li>
  </ul>
</p>
''').render(
    ATOM_CATEGORIES=ATOM_CATEGORIES.items(),
)

def test_canonical_rep():
    test_cases = (
        ('H,C4,H|SI|C|C2,H,C4', 'C4,H,H|SI|C|C4,C2,H'), # M(SI) > M(C)
        ('H,C4,H|C|C|C4,H', 'C4,H,H|C|C|C4,H'), #M(C) == M(C), 3 > 2
        ('H,C4,H|C|C|C3,H,H', 'C4,H,H|C|C|C3,H,H'), #M(C) == M(C), 3 ==3, 4 > 3
        ('O2|C|C|N2,H', 'N2,H|C|C|O2'), # M(O) > M(N)
        ('F,C4|C|C|O1,N3', 'F,C4|C|C|O1,N3'), # M(F) > M(O)
    )

    for (dihedral_string, solution) in test_cases:
        answer = str(FragmentDihedral(dihedral_string))
        assert answer == solution, '"{0}" (answer) != "{1}" (expected)'.format(answer, solution)

        cyclic_answer = str(FragmentDihedral(answer))
        assert cyclic_answer == answer, '"{0}" (answer) != "{1}" (expected)'.format(cyclic_answer, answer)

def test_patterns():
    test_cases = (
        'X,X|C|C|X,X',
        '!X+|N|C|CX',
        'C,X,B,D|N|N|H,_,C',
        'C+|CC|C|C+',
        'CL,CL,X{2}|C|Z|CX',
        'J+|C|S|H',
        '!J{2-4}|C|C|C',
        'X{2},I|C|Z|%',
        '!X{2},I|C|Z|%',
        'CL{2-3}|C|C|BR{3-5}',
    )

    for pattern in test_cases:
        print pattern
        print sql_pattern_matching_for(pattern)
        print


if __name__ == "__main__" :
    if True:
        test_canonical_rep()

    if False:
        test_patterns()

    if False:
        dihedral_1 = FragmentDihedral("C,C,H|C|C|C,H,H")
        print dihedral_1
        dihedral_2 = FragmentDihedral("C,H,H|C|C|C,C,H")
        print dihedral_2
        print dihedral_1 == dihedral_2
        dihedral_3 = FragmentDihedral(atom_list=(['H', 'H'], 'C', 'C', ['Cl', 'Cl']))
        print dihedral_3
        print dihedral_3 == dihedral_2

    assert re_pattern_matching_for('Z,%|Z|Z|Z,%', debug=True)('C,H|C|C|C,H') == True
    assert re_pattern_matching_for('Z|Z|Z|Z,%', debug=True)('C,H|C|C|C,H') == False
    assert re_pattern_matching_for('Z|Z|Z|Z,%', debug=True)('C|C|C|C,H') == True
    assert re_pattern_matching_for('N,J|C|C|J{3}', debug=True)('N,H|C|C|C,H,H') == True
    assert re_pattern_matching_for('N,J|C|C|J{2}', debug=True)('N,H|C|C|C,H') == True

    print re_pattern_matching_for('Z|Z|Z|Z', debug=True)

    print sql_pattern_matching_for('J{3}|C|C|J{3}')

