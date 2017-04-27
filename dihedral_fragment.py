from copy import deepcopy, copy
from itertools import product
from itertools import permutations
from jinja2 import Template
from re import search, sub, findall
from collections import namedtuple
from operator import itemgetter
from typing import Optional, Any, Tuple, Union, Sequence, NamedTuple, List, Callable

from atb_helpers.iterables import group_by
from atb_helpers.elements import ELEMENT_NUMBERS

Dihedral_Fragment_Str = str

Regular_Expression = str

Dihedral_Matching_Pattern = str

DEBUG = False

ENFORCE_CANONICAL_TRICYLIC = False

def print_if_DEBUG(something: Any) -> None:
    if DEBUG:
        print(something)

NEIGHBOUR_SEPARATOR = ','

def join_neighbours(neighbours: List[str]) -> str:
    return NEIGHBOUR_SEPARATOR.join(neighbours)

def split_neighbour_str(neighbour_str: str) -> List[str]:
    return neighbour_str.split(NEIGHBOUR_SEPARATOR)

GROUP_SEPARATOR = '|'

def join_groups(groups: List[str]) -> str:
    return GROUP_SEPARATOR.join(groups)

def split_group_str(group_str: str) -> List[str]:
    return group_str.split(GROUP_SEPARATOR)

class Dihedral_Fragment_Matching_Pattern(object):
    def __init__(self, pattern_string: Optional[str] = None) -> None:
        splitted_string = split_group_str(pattern_string)
        self.neighbours_1 = [atom.upper() for atom in split_neighbour_str(splitted_string[0])]
        self.atom_2 = splitted_string[1].upper()
        self.atom_3 = splitted_string[2].upper()
        self.neighbours_4 = [atom.upper() for atom in split_neighbour_str(splitted_string[3])]

def has_regex_pattern(pattern: str) -> bool:
    assert ('$' not in pattern) and ('^' not in pattern), pattern
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

def REGEX_NOT(pattern: str) -> str:
    return '[^(' + str(pattern) + ')]'

def REGEX_GROUP(index: int) -> str:
    return BACKSLASH + str(index)

def REGEX_OR(*args: List[str]) -> str:
    return '(' + REGEX_OR_OPERATOR.join(args) + ')'

def REGEX_ESCAPE(pattern: str, flavour: str = 'sql') -> str:
    if flavour == 'sql':
        return BACKSLASH + str(pattern)
        return BACKSLASH + BACKSLASH + str(pattern)
    elif flavour == 're':
        return '[' + str(pattern) + ']'
    else:
        raise Exception()

def REGEX_AT_LEAST(pattern: str, escape_plus: bool = True) -> str:
    return str(pattern) + (BACKSLASH if escape_plus else '') + '+'

def FORMAT_ESCAPED(pattern: str) -> str:
    return pattern.replace('{', '{{').replace('}', '}}')

def FORMAT_UNESCAPED(pattern: str) -> str:
    return pattern.replace('{{', '{').replace('}}', '}')

def NOT(pattern: str) -> str:
    return '!' + str(pattern)

def CAPTURE(pattern: str) -> str:
    return '(' + str(pattern) + ')'

def ESCAPE(pattern: str) -> str:
    return BACKSLASH + str(pattern)

NO_VALENCE = None

def element_valence_for_atom(atom_desc: str) -> Tuple[str, Optional[int]]:
    upper_atom = atom_desc.upper()
    match = search(
        CAPTURE('[' + ATOM_CHARACTERS + ']') + CAPTURE('[' + VALENCE_CHARACTERS + ']'),
        upper_atom,
    )
    if match:
        element, valence_str = match.groups()
        valence = int(valence_str)
    else:
        element, valence = upper_atom, NO_VALENCE
    return (element, valence)

def on_asc_number_electron_then_asc_valence(atom) -> Union[Tuple[int, int], Tuple[int]]:
    ASC = lambda x: x
    element, valence = element_valence_for_atom(atom)
    try:
        return (
            ASC(ELEMENT_NUMBERS[element]),
            ASC(valence),
        )
    except KeyError:
        return (
            999,
        )

Cycle = NamedTuple('Cycle', [('i', int), ('n', int), ('j', int)])

def Small_Cycle(*args: Sequence[Any]) -> Cycle:
    assert len(args) == 3, 'Cycles of length > 9 are not allowed.'
    return Cycle(*args)

GROUP_INDICES = (0, 1, 2, 3, 4)
LEFT_GROUP_INDEX, LEFT_ATOM_INDEX, RIGHT_ATOM_INDEX, RIGHT_GROUP_INDEX, CYCLES_INDEX = GROUP_INDICES

CHIRAL_MARKER = '*'

class Dihedral_Fragment(object):
    def __init__(self, dihedral_string: Optional[str] = None, atom_list: Optional[Any] = None) -> None:
        assert dihedral_string is not None or atom_list is not None, [dihedral_string, atom_list]

        if dihedral_string is not None:
            splitted_string = split_group_str(dihedral_string)
            self.neighbours_1 = [atom.upper() for atom in split_neighbour_str(splitted_string[LEFT_GROUP_INDEX])]
            self.atom_2 = splitted_string[LEFT_ATOM_INDEX].upper()
            self.atom_3 = splitted_string[RIGHT_ATOM_INDEX].upper()
            self.neighbours_4 = [atom.upper() for atom in split_neighbour_str(splitted_string[RIGHT_GROUP_INDEX])]
            self.cycles = (
                [
                    Small_Cycle(*list(map(int, cycle_str)))
                    for cycle_str in
                    split_neighbour_str(splitted_string[CYCLES_INDEX])
                ]
                if len(splitted_string) == 5
                else []
            )
        else:
            if len(atom_list) == 4:
                self.neighbours_1, self.atom_2, self.atom_3, self.neighbours_4 = atom_list
                self.cycles = []
            elif len(atom_list) == 5:
                self.neighbours_1, self.atom_2, self.atom_3, self.neighbours_4, self.cycles = atom_list
                self.cycles = [Small_Cycle(*c) for c in self.cycles]
            else:
                raise Exception('Wrong length of atom_list: {0}'.format(atom_list))

        canonical_rep = self.__canonical_rep__()

        self.neighbours_1 = canonical_rep.neighbours_1
        self.atom_2 = canonical_rep.atom_2
        self.atom_3 = canonical_rep.atom_3
        self.neighbours_4 = canonical_rep.neighbours_4
        self.cycles = canonical_rep.cycles

    def has_cycles(self) -> bool:
        return bool(self.cycles)

    def __str__(self, flag_chiral_sides: bool = False) -> str:
        return join_groups(
            [
                join_neighbours(self.neighbours_1),
                self.atom_2
                +
                ('' if (not flag_chiral_sides or (flag_chiral_sides and not self.is_left_chiral())) else CHIRAL_MARKER),
                self.atom_3
                +
                ('' if (not flag_chiral_sides or (flag_chiral_sides and not self.is_right_chiral())) else CHIRAL_MARKER),
                join_neighbours(self.neighbours_4),
            ]
            +
            (
                [
                    ','.join(
                        list(map(
                            lambda cycle: ''.join(map(str, cycle)),
                            self.cycles,
                        )),
                    ),
                ]
                if self.has_cycles()
                else []
            )
        )

    def __eq__(self, other: Any) -> bool:
        return self.__dict__ == other.__dict__
        #return self.__canonical_rep__().__dict__ == other.__canonical_rep__().__dict__

    def __ne__(self, other: Any) -> bool:
        return not self == other

    def __canonical_rep__(self) -> Any:
        other = copy(self)

        def sorted_neighbours_permutation_dict(neighbours):
            sorted_neighbours = list(
                sorted(
                    enumerate(neighbours),
                    key=lambda i_neighbour: on_asc_number_electron_then_asc_valence(i_neighbour[1]),
                    reverse=True,
                )
            )
            permutation_dict = dict(
                [(i, j) for (j, (i, _)) in enumerate(sorted_neighbours)]
            )
            return (
                list(map(itemgetter(1), sorted_neighbours)),
                permutation_dict,
            )


        def sort_neighbours_renumber_cycles(fragment):
            '''Order each neighbour list by alphabetical order, reordering the (possible) cycles to match those changes.'''
            # WARNING: Keep in mind that maintaining order will be necessary for maintaining sterochemistry information
            fragment.neighbours_1, permutation_1 = sorted_neighbours_permutation_dict(fragment.neighbours_1)
            fragment.neighbours_4, permutation_4 = sorted_neighbours_permutation_dict(fragment.neighbours_4)
            fragment.cycles = [
                Cycle(permutation_1[neighbour_id_1], cycle_length, permutation_4[neighbour_id_4])
                for (neighbour_id_1, cycle_length, neighbour_id_4) in self.cycles
            ]

        def canonise_cycles(fragment):
            if len(fragment.cycles) in (0, 1):
                pass
            elif len(fragment.cycles) == 2:
                cycle_0, cycle_1 = fragment.cycles
                should_order_left = (fragment.neighbours_1[cycle_0.i] == fragment.neighbours_1[cycle_1.i])
                should_order_right = (fragment.neighbours_4[cycle_0.j] == fragment.neighbours_4[cycle_1.j])
                fragment.cycles = [
                    Cycle(
                        min(cycle_0.i, cycle_1.i) if should_order_left else cycle_0.i,
                        min(cycle_0.n, cycle_1.n),
                        min(cycle_0.j, cycle_1.j) if should_order_right else cycle_0.j,
                    ),
                    Cycle(
                        max(cycle_0.i, cycle_1.i) if should_order_left else cycle_1.i,
                        max(cycle_0.n, cycle_1.n),
                        max(cycle_0.j, cycle_1.j) if should_order_right else cycle_1.j,
                    ),
                ]
            else:
                Is, Ns, Js = list(zip(*self.cycles))
                should_order = dict(
                    [
                        (direction, len(set(neighbours)) == 1)
                        for (direction, neighbours) in zip(
                            ('left', 'right'),
                            (fragment.neighbours_1, fragment.neighbours_4),
                        )
                    ],
                )

                if ENFORCE_CANONICAL_TRICYLIC:
                    assert all(should_order.values()), 'Only symmetric (all similar or all different) environments are allowed for N-cycles where N >= 3 (fragment={0}, cycles={1}).'.format(
                        str(fragment),
                        fragment.cycles,
                    )

                fragment.cycles = [
                    Cycle(i, n, j)
                    for (i, n, j) in
                    zip(
                        *list(map(
                            sorted,
                            (Is, Ns, Js),
                        ))
                    )
                ]

                def cycle_key(cycle):
                    return(
                        fragment.neighbours_1[cycle.i],
                        len([1 for other_cycle in fragment.cycles if cycle.i == other_cycle.i]),
                        fragment.neighbours_4[cycle.j],
                        len([1 for other_cycle in fragment.cycles if cycle.j == other_cycle.j]),
                    )

        def order_cycles(fragment):
            fragment.cycles.sort(
                key=lambda cycle: (cycle.n, cycle.i, cycle.j),
            )

        sort_neighbours_renumber_cycles(other)
        canonise_cycles(other)
        order_cycles(other)

        # Compare the two central atoms and put the heavier one on the left
        if on_asc_number_electron_then_asc_valence(other.atom_2) <  on_asc_number_electron_then_asc_valence(other.atom_3):
            should_reverse = True
        elif on_asc_number_electron_then_asc_valence(other.atom_2) == on_asc_number_electron_then_asc_valence(other.atom_3):
            should_reverse = False

            if len(other.neighbours_1) < len(other.neighbours_4):
                should_reverse = True
            elif len(other.neighbours_1) == len(other.neighbours_4):
                # If identical central atoms, and same number of neighbours on both ends, try to resolve ambiguity one neighbour at a time
                for (neighbour_1, neighbour_4) in zip(other.neighbours_1, other.neighbours_4):
                    if on_asc_number_electron_then_asc_valence(neighbour_1) < on_asc_number_electron_then_asc_valence(neighbour_4):
                        should_reverse = True
                        break
                    elif on_asc_number_electron_then_asc_valence(neighbour_1) == on_asc_number_electron_then_asc_valence(neighbour_4):
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

    def reverse_dihedral(self) -> None:
        self.atom_2, self.atom_3 = self.atom_3, self.atom_2
        self.neighbours_1, self.neighbours_4 = self.neighbours_4, self.neighbours_1
        self.cycles = [x[::-1] for x in self.cycles]

    def is_left_chiral(self) -> bool:
        return len(set(self.neighbours_1)) == 3

    def is_right_chiral(self) -> bool:
        return len(set(self.neighbours_4)) == 3

    def is_chiral_fragment(self) -> bool:
        return (self.is_left_chiral() or self.is_right_chiral())

Operator_Pattern = NamedTuple('Operator_Pattern', [('pattern', str), ('replacement', str), ('substitution_type', str)])

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

SQL_SUBSTITUTION_CHARACTERS = ('_', '%')
SQL_FULL_REGEX_CHARACTERS = ('{', '}', '!', '+', '-') + tuple(ATOM_CATEGORIES.keys())
has_substitution_pattern = lambda x: any([y in x for y in SQL_SUBSTITUTION_CHARACTERS])

def substitute_atom_pattern(atom_pattern):
    def cleaned_atom_pattern(x: str) -> str:
        return sub('[\[\]\(\)\^\|\{\}\+\-0-9,]', '', x) # This line prevents atoms from ever having numbers in their name

    for cleaned_atom_pattern in (atom_pattern, cleaned_atom_pattern(atom_pattern)):
        if cleaned_atom_pattern in ATOM_CATEGORIES:
            return sub(cleaned_atom_pattern, ATOM_CATEGORIES[cleaned_atom_pattern], atom_pattern)
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

def split_on_atoms(group: Any) -> Any:
    atom_patterns = split_neighbour_str(group)
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

    left_neighbour_groups, right_neighbour_groups = list(map(
        lambda index: split_neighbour_str(components[index]),
        (LEFT_GROUP_INDEX, RIGHT_GROUP_INDEX),
    ))

    left_permutations, right_permutations = list(map(
        lambda group: permutations(list(range(len(group_by(group, key=lambda x:x))))),
        (left_neighbour_groups, right_neighbour_groups),
    ))

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
        list(zip(
            sorted(group_by(component_list, lambda x:x).keys()),
            permutation,
        ))
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

def test_canonical_rep() -> None:
    test_cases = (
        ('H,C4,H|SI|C|C2,H,C4', 'C4,H,H|SI|C|C4,C2,H'), # M(SI) > M(C)
        ('H,C4,H|C|C|C4,H', 'C4,H,H|C|C|C4,H'), #M(C) == M(C), 3 > 2
        ('H,C4,H|C|C|C3,H,H', 'C4,H,H|C|C|C3,H,H'), #M(C) == M(C), 3 ==3, 4 > 3
        ('O2|C|C|N2,H', 'N2,H|C|C|O2'), # M(O) > M(N)
        ('F,C4|C|C|O1,N3', 'F,C4|C|C|O1,N3'), # M(F) > M(O)
    )

    for (dihedral_string, solution) in test_cases:
        answer = str(Dihedral_Fragment(dihedral_string))
        assert answer == solution, '"{0}" (answer) != "{1}" (expected)'.format(answer, solution)

        cyclic_answer = str(Dihedral_Fragment(answer))
        assert cyclic_answer == answer, '"{0}" (answer) != "{1}" (expected)'.format(cyclic_answer, answer)

def test_patterns() -> None:
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
        print(pattern)
        print(sql_pattern_matching_for(pattern))
        print()

def test_cyclic_fragments() -> None:
    cyclic_fragment = Dihedral_Fragment(atom_list=(['H', 'H', 'C'], 'C', 'C', ['H', 'C', 'H'], [[2, 0, 1]]))
    print(str(cyclic_fragment))
    assert str(cyclic_fragment) == 'C,H,H|C|C|C,H,H|000', cyclic_fragment

    polycyclic_fragment, answer = Dihedral_Fragment(atom_list=(['H', 'C', 'C'], 'C', 'C', ['C', 'C', 'H'], [[2, 2, 1], [1, 3, 0]])), 'C,C,H|C|C|C,C,H|020,131'
    print(str(polycyclic_fragment))
    assert str(polycyclic_fragment) == answer, '{0} != {1} (expected)'.format(
        str(polycyclic_fragment),
        answer,
    )
    print(Dihedral_Fragment(str(polycyclic_fragment)))

    polycyclic_fragment, answer = Dihedral_Fragment(atom_list=(['H', 'N', 'O'], 'C', 'C', ['C', 'C', 'H'], [[2, 2, 1], [1, 3, 0]])), 'C,C,H|C|C|O,N,H|020,131'
    print(str(polycyclic_fragment))
    assert str(polycyclic_fragment) == answer, '{0} != {1} (expected)'.format(
        str(polycyclic_fragment),
        answer,
    )
    print(Dihedral_Fragment(str(polycyclic_fragment)))

    assert str(Dihedral_Fragment('C,C,C|C|C|C,C,C|002,101,200')) == 'C,C,C|C|C|C,C,C|000,101,202'

    try:
        Dihedral_Fragment('C,C,N|C|C|C,C,C|002,101,200')
        raise Exception('This should have failed.')
    except AssertionError:
        print('Non-symmetric N=3 rings failed as expected.')

    try:
        Dihedral_Fragment('C,C,C|C|C|C,C,C|0100')
        raise Exception('This should have failed.')
    except AssertionError:
        print('Rings length > 9 failed as expected.')

def test_atom_list_init() -> None:
    fragment = Dihedral_Fragment(atom_list=(['C'], 'C', 'C', ['C']))

def test_misc() -> None:
    dihedral_1 = Dihedral_Fragment("C,C,H|C|C|C,H,H")
    print(dihedral_1)
    dihedral_2 = Dihedral_Fragment("C,H,H|C|C|C,C,H")
    print(dihedral_2)
    print(dihedral_1 == dihedral_2)
    dihedral_3 = Dihedral_Fragment(atom_list=(['H', 'H'], 'C', 'C', ['Cl', 'Cl']))
    print(dihedral_3)
    print(dihedral_3 == dihedral_2)

def test_chiral_str() -> None:
    dihedral_1 = Dihedral_Fragment("C,C,H|C|C|C,H,H")
    print(dihedral_1.__str__())
    print(dihedral_1.__str__(flag_chiral_sides=True))

if __name__ == "__main__" :
    #test_atom_list_init()
    test_canonical_rep()
    test_chiral_str()
    test_cyclic_fragments()
    test_misc()

    if False:
        test_patterns()

    assert re_pattern_matching_for('Z,%|Z|Z|Z,%', debug=True)('C,H|C|C|C,H') == True
    assert re_pattern_matching_for('Z|Z|Z|Z,%', debug=True)('C,H|C|C|C,H') == False
    assert re_pattern_matching_for('Z|Z|Z|Z,%', debug=True)('C|C|C|C,H') == True
    assert re_pattern_matching_for('N,J|C|C|J{3}', debug=True)('N,H|C|C|C,H,H') == True
    assert re_pattern_matching_for('N,J|C|C|J{2}', debug=True)('N,H|C|C|C,H') == True

    print(re_pattern_matching_for('Z|Z|Z|Z', debug=True))

    print(sql_pattern_matching_for('J{3}|C|C|J{3}'))

