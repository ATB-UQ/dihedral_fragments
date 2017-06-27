from re import search, sub, findall
from operator import itemgetter
from typing import List, Tuple, Sequence, Dict, Callable, Any, NamedTuple, Optional
from itertools import product, permutations, groupby
from jinja2 import Template

from dihedral_fragments.dihedral_fragment import Dihedral_Fragment, split_neighbour_str, split_group_str, LEFT_ATOM_INDEX, RIGHT_ATOM_INDEX, LEFT_GROUP_INDEX, RIGHT_GROUP_INDEX, join_groups, join_neighbours, print_if_DEBUG, DEBUG, Dihedral_Fragment_Str
from dihedral_fragments.regex import CAPTURE, NOT, ESCAPE, exactly_N_times_operator, N_to_M_times_operator, REGEX_START_ANCHOR, REGEX_END_ANCHOR, REGEX_NOT_SET, REGEX_GROUP, REGEX_SET, ESCAPED_COMMA, UNESCAPE_COMMA, ONE_ATOM, REGEX_AT_LEAST, REGEX_OR, ANY_NUMBER_OF_ATOMS, REGEX_ESCAPE, FORMAT_ESCAPED, FORMAT_UNESCAPED

Operator_Pattern = NamedTuple('Operator_Pattern', [('pattern', str), ('replacement', str), ('substitution_type', str)])

Dihedral_Matching_Pattern = str

Regular_Expression = str

class Dihedral_Fragment_Matching_Pattern(object):
    def __init__(self, pattern_string: Optional[str] = None) -> None:
        splitted_string = split_group_str(pattern_string)
        self.neighbours_1 = [atom.upper() for atom in split_neighbour_str(splitted_string[0])]
        self.atom_2 = splitted_string[1].upper()
        self.atom_3 = splitted_string[2].upper()
        self.neighbours_4 = [atom.upper() for atom in split_neighbour_str(splitted_string[3])]

OPERATORS = (
    (   # !A
        NOT(CAPTURE(ONE_ATOM)),
        REGEX_NOT_SET( REGEX_GROUP(1) ),
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

OPERATOR_PATTERNS = list(map(
    lambda x: Operator_Pattern(*x),
    OPERATORS,
))

def regex_filters(flavour: str ='sql') -> List[Any]:
    return [
        # Order matters here !
        Operator_Pattern(
            '|',
            REGEX_ESCAPE('|', flavour=flavour),
            'str',
        ),
    ] + OPERATOR_PATTERNS + [
        Operator_Pattern('%', ANY_NUMBER_OF_ATOMS, 'str'),
    ]

ATOM_CATEGORIES = {
    'J': REGEX_OR('C', 'H'),
    'X': REGEX_OR('F', 'I', 'BR', 'CL'),
    'Y': REGEX_OR('N', 'O', 'S'),
    'Z': FORMAT_ESCAPED(ONE_ATOM),
}

SQL_SUBSTITUTION_CHARACTERS = ['_', '%']
SQL_FULL_REGEX_CHARACTERS = ['{', '}', '!', '+', '-'] + list(ATOM_CATEGORIES.keys())
has_substitution_pattern = lambda x: any([y in x for y in SQL_SUBSTITUTION_CHARACTERS])

def has_regex_pattern(pattern: str) -> bool:
    assert ('$' not in pattern) and ('^' not in pattern), pattern
    return any([y in pattern for y in SQL_FULL_REGEX_CHARACTERS])

def substitute_atom_pattern(atom_pattern: str) -> str:
    def cleaned_atom_pattern(x: str) -> str:
        '''This line prevents atoms from ever having numbers in their name'''
        return sub(
            REGEX_SET(r'\\', '\[', '\]', '\(', '\)', '\^', '\|', '\{', '\}', '\+', '\-', '0-9', ','),
            '',
            x,
        )

    all_atom_patterns = [atom_pattern, cleaned_atom_pattern(atom_pattern)]

    if DEBUG:
        print('all_atom_patterns', all_atom_patterns)

    for an_atom_pattern in all_atom_patterns:
        if an_atom_pattern in ATOM_CATEGORIES:
            return sub(an_atom_pattern, ATOM_CATEGORIES[an_atom_pattern], atom_pattern)
    else:
        return atom_pattern

def substitute_atoms_in_pattern(pattern: str, flavour: str = 'sql') -> str:
    if flavour == 'sql':
        BAR = REGEX_ESCAPE('|')
    elif flavour == 're':
        BAR = '|'
    return BAR.join(
        [
            join_neighbours(
                [
                    substitute_atom_pattern(atom)
                    for atom in split_on_atoms(group)
                ],
            )
            for group in split_group_str(pattern)
        ]
    )

def split_on_atoms(group: Any) -> Any:
    atom_patterns = split_neighbour_str(group)
    if DEBUG:
        print('group', group)
        print('atom_patterns', atom_patterns)
        print()
    return atom_patterns

def apply_regex_filters(string, debug=False, flavour='sql'):
    for (pattern, replacement, substitution_type) in regex_filters(flavour):
        if debug:
            old_string = string

        if substitution_type == 'str':
            string = string.replace(pattern, replacement)
        elif substitution_type == 're':
            string = sub(pattern, replacement, string)
        elif substitution_type == 're_map_groups':
            general_pattern = pattern()

            matches = findall(general_pattern, string)
            if matches:
                if debug:
                    print(matches)
                for match in matches:
                    tailored_pattern = pattern(*match)
                    mapped_replacement = replacement(*match)
                    if debug:
                        print(tailored_pattern, match)
                        print(mapped_replacement)
                    string = sub(tailored_pattern, mapped_replacement, string)
        else:
            raise Exception('Wrong substitution_type')

        if debug:
            print([pattern, replacement, old_string], string)
    print_if_DEBUG(string)
    return UNESCAPE_COMMA(substitute_atoms_in_pattern(string, flavour=flavour))

def escaped_special_regex_characters(patterns, flavour='sql'):
    return [ (REGEX_START_ANCHOR + apply_regex_filters(pattern, flavour=flavour)  + REGEX_END_ANCHOR) for pattern in patterns]

def sql_OR(*args: Sequence[str]) -> str:
    return ' '.join(
        ['('] + [' OR '.join(*args)] + [')']
    )

TRUE_OR_FALSE = [False, True]
FALSE = [False]
PUT_SUBSTITUTION_PATTERN_FIRST = True


ON_SELF = lambda x: x

def re_patterns(pattern: Dihedral_Matching_Pattern, full_regex: bool = False, flavour: str = 'sql', debug: bool = False, metadata: Any = None) -> List[Regular_Expression]:
    components = split_group_str(pattern)
    assert len(components) == 4
    need_to_reverse_inner_atoms = (components[LEFT_ATOM_INDEX] == components[RIGHT_ATOM_INDEX]) or any([x in list(ATOM_CATEGORIES.keys()) for x in (components[LEFT_ATOM_INDEX], components[RIGHT_ATOM_INDEX])])

    if debug:
        print()
        print(metadata)
        print('need_to_reverse_inner_atoms: {0}'.format(need_to_reverse_inner_atoms))
        print(components[LEFT_ATOM_INDEX])
        print(components[RIGHT_ATOM_INDEX])

    left_neighbour_groups, right_neighbour_groups = list(
        map(
            lambda index: split_neighbour_str(components[index]),
            (LEFT_GROUP_INDEX, RIGHT_GROUP_INDEX),
        ),
    )

    left_permutations, right_permutations = list(
        map(
            lambda group: permutations(
                range(
                    len(
                        list(groupby(sorted(group, key=ON_SELF), key=ON_SELF)),
                    ),
                ),
            ),
            (left_neighbour_groups, right_neighbour_groups),
        ),
    )

    patterns = [
        correct_pattern(
            pattern,
            should_reverse,
            left_permutation,
            right_permutation,
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

def re_pattern_matching_for(pattern: Dihedral_Matching_Pattern, debug: bool = False, metadata: Any = None) -> Callable[[Dihedral_Fragment_Str], bool]:
    if not (has_substitution_pattern(pattern) or has_regex_pattern(pattern)):
        return lambda test_string: test_string == str(Dihedral_Fragment(pattern))
    else:
        patterns = [FORMAT_UNESCAPED(re_pattern) for re_pattern in re_patterns(pattern, full_regex=True, flavour='re', debug=debug, metadata=metadata)]

        def match_pattern_to(test_string):
            if debug:
                print(metadata)
                print(pattern)
                print(patterns)
                print(test_string)
                print(any([search(match_pattern, test_string) for match_pattern in patterns]))
                print()
            return any([search(match_pattern, test_string) for match_pattern in patterns])

        return match_pattern_to

def sql_pattern_matching_for(pattern: Dihedral_Matching_Pattern, matching_field_name: str = 'dihedral_string'):
    if not (has_substitution_pattern(pattern) or has_regex_pattern(pattern)):
        return '{matching_field_name}="{pattern}"'.format(
            matching_field_name=matching_field_name,
            pattern=Dihedral_Fragment(pattern),
        )
    else:
        use_full_regex = has_regex_pattern(pattern)
        all_patterns = re_patterns(pattern, full_regex=use_full_regex, flavour='sql')
        #patterns = [FORMAT_UNESCAPED(re_pattern) for re_pattern in re_patterns(pattern, full_regex=True, flavour='sql')]

        return sql_OR(
            [
                '{matching_field_name} {sql_operator} "{match_pattern}"'.format(
                    matching_field_name=matching_field_name,
                    sql_operator=('LIKE' if not use_full_regex else 'REGEXP'),
                    match_pattern=match_pattern,
                )
                for match_pattern in all_patterns
            ]
        )

def correct_pattern(pattern: Dihedral_Matching_Pattern, should_reverse: bool, left_permutation: Sequence[int], right_permutation: Sequence[int]) -> str:
    return join_groups(
        (reversed if should_reverse else (lambda x: x))(
            [
                sorted_components(
                    i,
                    x,
                    left_permutation,
                    right_permutation,
                )
                for (i, x) in enumerate(split_group_str(pattern))
            ]
        )
    )

def sorted_components(component_index: int, component: str, left_permutation: Sequence[int], right_permutation: Sequence[int]) -> str:
    if component_index in (0,3):
        return join_neighbours(
            sorted_components_list(
                split_neighbour_str(component),
                permutation=left_permutation if component_index == 0 else right_permutation,
            )
        )
    else:
        return component

def sorted_components_list(component_list: List[str], permutation: Sequence[int]) -> List[str]:
    sorting_dict = dict(
        zip(
            sorted(
                [
                    key
                    for (key, group) in groupby(
                        sorted(component_list, key=ON_SELF),
                        key=ON_SELF,
                    )
                ],
            ),
            permutation,
        ),
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
      <li><code>{{ code }}</code>Any (single) atom in {{ in_code_tag(matches) }}</li>
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
    ATOM_CATEGORIES=list(ATOM_CATEGORIES.items()),
    in_code_tag=lambda x: '<code>' + (str(x) if x else '') + '</code>',
)


