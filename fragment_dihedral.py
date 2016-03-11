from copy import deepcopy, copy
from itertools import product
from itertools import permutations
from jinja2 import Template
import re

from atb_helpers.iterables import group_by
from pipeline.pipelineHelperFunctions import ELEMENT_NUMBERS

CHEMICAL_GROUPS = (
    # Hydrocarbons
    ('alkane', 'J,J,J|C|C|J,J,J'),
    ('alkane', 'J,J|C|C|J,J'),
    ('alkyne', 'J|C|C|J'),

    # Groups containing halogens
    ('monofluoro', 'F,!X,!X|C|Z|%'),
    ('difluoro', 'F,F,!X|C|Z|%'),
    ('trifluoro', 'F,F,F|C|Z|%'),
    ('chloro', 'CL,!X,!X|C|Z|%'),
    ('dichloro', 'CL,CL,!X|C|Z|%'),
    ('trichloro', 'CL,CL,CL|C|Z|%'),
    ('bromo', 'BR,!X,!X|C|Z|%'),
    ('dibromo', 'BR,BR,!X|C|Z|%'),
    ('tribromo', 'BR,BR,BR|C|Z|%'),
    ('iodo', 'I,!X,!X|C|Z|%'),
    ('diiodo', 'I,I,!X|C|Z|%'),
    ('triiodo', 'I,I,I|C|Z|%'),

    # Groups containing oxygen
    ('alcohol I', 'C,H,H|C|O|H'),
    ('alcohol II', 'C,C,H|C|O|H'),
    ('alcohol III', 'C,C,C|C|O|H'),
    ('ketone', '%|C|C|O,C'),
    ('aldehyde', '%|C|C|O,H'),
    ('acyl halide', '%|C|C|O,X'),
    ('carboxylic acid', '%,O|C|O|H'),
    ('ester', '%,O|C|O|C'),

    # Groups containing nitrogen
    ('carboxamide', ''),
    ('amine I', '%|C|N|H,H'),
    ('amine II', '%|C|N|C,H'),
    ('amine III', '%|C|N|C,C'),
    ('ammonium ion', '%|J|N|J,J,J'),
    ('ketimine I', 'H|N|C|C,C'),
    ('ketimine II', 'C|N|C|C,C'),
    ('aldimine I', 'H|N|C|C,H'),
    ('aldimine II', 'C|N|C|C,H'),
    ('imide', 'N/A'),
    ('azide', 'N|N|N|C'),
    ('azo', 'C|N|N|C'),
    ('cyanate', 'C|O|C|N'),
    ('isocyanate', 'C|N|C|O'),
    ('nitrate', 'C|O|N|O,O'),
    ('nitrile', 'N|C|C|Z,Z,Z'),
    ('isonitrile', 'C|N|C|Z,Z,Z'),
    ('nitro', 'C|N|O,O'),
    ('nitroso', '%|C|N|O'),
    ('pyridyl', ''),

    # Groups containing sulfur
    ('thiol', 'J+|C|S|H'),
    ('sulfide', 'J+|C|S|C'),
    ('disulfide', 'C|S|S|C'),
    ('sulfinyl', '%|C|S|C,O'),
    ('sulfonyl', '%|C|S|C,O,O'),
    ('sulfino', 'C,O|S|O|H'),
    ('sulfo', 'C,O,O|S|O|H'),
    ('thiocyanate', 'C|S|C|N'),
    ('isothiocyanate', 'C|N|C|S'),
    ('carbonothioyl', '%|J|C|S,J'),

    # Groups containing phosphorus
    ('phosphino', '%|C|P|C,C'),
    ('phosphono', '%|C|P|O,O,O'),
    ('phosphate', 'C|O|P|O,O,O'),

    # Groups containing boron
    ('borono', '%|C|B|O,O'),
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

def has_regex_pattern(pattern):
    assert ('$' not in pattern) and ('^' not in pattern)
    return any([y in pattern for y in SQL_FULL_REGEX_CHARACTERS])

REGEX_START, REGEX_END = ('^', '$')

BACKSLASH = '\\'
COMMA = ','
ONE_LETTER = '[A-Z]'
ONE_ATOM = ONE_LETTER + '{1,2}'
ONE_NUMBER = '[0-9]'
ANY_NUMBER_OF_ATOMS = '[A-Z,]*'

def REGEX_NOT(pattern):
    return '[^(' + str(pattern) + ')]'

def REGEX_GROUP(index):
    return BACKSLASH + str(index)

def REGEX_OR(*args):
    return '(' + '|'.join(args) + ')'

def REGEX_ESCAPE(pattern):
    return BACKSLASH + BACKSLASH + str(pattern)

def REGEX_AT_LEAST(pattern, escape_plus=True):
    return str(pattern) + (BACKSLASH if escape_plus else '') + '+'

def FORMAT_ESCAPED(pattern):
    return pattern.replace('{', '{{').replace('}', '}}')

def NOT(pattern):
    return '!' + str(pattern)

def CAPTURE(pattern):
    return '(' + str(pattern) + ')'

def ESCAPE(pattern):
    return BACKSLASH + str(pattern)

OPERATORS = (
    (   # !A
        NOT(CAPTURE(ONE_ATOM)),
        REGEX_NOT( REGEX_GROUP(1) ),
        're',
    ),
    (   # A+
        REGEX_AT_LEAST(CAPTURE(ONE_ATOM), escape_plus=True),
        REGEX_AT_LEAST(REGEX_OR( REGEX_GROUP(1), COMMA), escape_plus=False),
        're',
    ),
    (   # A{X}
        CAPTURE(ONE_ATOM) + ESCAPE('{') + CAPTURE(ONE_NUMBER) + ESCAPE('}'),
        lambda x, z: FORMAT_ESCAPED( REGEX_OR(x, COMMA) + '{' + str(2*int(z) - 1) + '}' ),
        're_map_groups',
    ),
    (   # A{X,Y}
        CAPTURE(ONE_ATOM) + ESCAPE('{') + CAPTURE(ONE_NUMBER) + ESCAPE('-') + CAPTURE(ONE_NUMBER) + ESCAPE('}'),
        lambda x, y, z: FORMAT_ESCAPED( REGEX_OR(x, COMMA) + '{' + str(int(y) + 1) + ',' + str(2*int(z) - 1) + '}' ),
        're_map_groups',
    ),
)

REGEX_FILTERS = (
    # Order matters here !
    (
        '|',
        REGEX_ESCAPE('|'),
        'str',
    ),
) + OPERATORS + (
    ('%', ANY_NUMBER_OF_ATOMS, 'str'),
)

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
        return re.sub('[\[\]\(\)\^\|\{\}\+\-0-9,]', '', x) # This line prevents atoms for ever having numbers in their name

    for cleaned_atom_pattern in (atom_pattern, cleaned_atom_pattern(atom_pattern)):
        if cleaned_atom_pattern in ATOM_CATEGORIES:
            return re.sub(cleaned_atom_pattern, ATOM_CATEGORIES[cleaned_atom_pattern], atom_pattern)
    return atom_pattern

def substitute_atoms_in_pattern(pattern):
    BAR = REGEX_ESCAPE('|')
    return BAR.join(
        [
            ','.join([ substitute_atom_pattern(atom) for atom in split_on_atoms(group)])
            for group in pattern.split(BAR)
        ]
    )

def split_on_atoms(groups):
    atom_patterns = groups.split(',')

    if len(atom_patterns):
        return atom_patterns

    atoms = []
    acc = []

    def discharge_acc():
        if acc is not []:
            atoms.append(','.join(acc))

    for atom_pattern in atom_patterns:
        if not (re.search(ONE_LETTER, atom_pattern[0]) or re.search(ONE_LETTER, atom_pattern[-1])):
            acc.append(atom_pattern)
        else:
            discharge_acc()
            acc = []
            atoms.append(atom_pattern)
    discharge_acc()

    return atoms

def apply_regex_filters(string, debug=False):
    for (pattern, replacement, substitution_type) in REGEX_FILTERS:
        if debug:
            old_string = string

        if substitution_type == 'str':
            string = string.replace(pattern, replacement)
        elif substitution_type == 're':
            string = re.sub(pattern, replacement, string)
        elif substitution_type == 're_map_groups':
            matches = re.findall(pattern, string)
            if matches:
                if len(matches) > 1:
                    raise Exception('This feature is not supported yet.')
                mapped_replacement = replacement(*matches[0])
                string = re.sub(pattern, mapped_replacement, string)
        else:
            raise Exception('Wrong substitution_type')

        if debug:
            print [pattern, replacement, old_string], string
    return substitute_atoms_in_pattern(string)

def escaped_special_regex_characters(patterns):
    return [ (REGEX_START + apply_regex_filters(pattern)  + REGEX_END) for pattern in patterns]

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
        return 'dihedral_string="{pattern}"'.format(pattern=FragmentDihedral(pattern))
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

if __name__ == "__main__" :
    for pattern in ('X,X|C|C|X,X', '!X+|N|C|CX', 'C,X,B,D|N|N|H,_,C', 'C+|CC|C|C+', 'CL,CL,X{2}|C|Z|CX', 'J+|C|S|H', '!J{2-4}|C|C|C', 'X{2},I|C|Z|%', '!X{2},I|C|Z|%', 'CL{2-3}|C|C|BR{3-5}'):
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


