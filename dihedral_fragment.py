from copy import deepcopy, copy
from re import search
from typing import Optional, Any, Tuple, Union, Sequence, NamedTuple, List, Callable, Dict
from sys import stderr

from dihedral_fragments.deque import deque, Deque, rotated_deque, reversed_deque
from dihedral_fragments.atomic_numbers import ATOMIC_NUMBERS
from dihedral_fragments.regex import CAPTURE, ATOM_CHARACTERS, VALENCE_CHARACTERS

Dihedral_Fragment_Str = str

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

def ASC(x: Optional[int]) -> Optional[int]:
    if x is None:
        return x
    else:
        return x

def DESC(x: Optional[int]) -> Optional[int]:
    if x is None:
        return x
    else:
        return -x

def on_asc_atomic_number_then_asc_valence(atom_desc: str) -> Tuple[int, int]:
    element, valence = element_valence_for_atom(atom_desc)
    try:
        return (
            ASC(ATOMIC_NUMBERS[element]),
            ASC(valence),
        )
    except KeyError:
        return (
            999,
            999,
        )

def on_desc_atomic_number_then_desc_valence(atom_desc: str) -> Tuple[int, int]:
    element, valence = element_valence_for_atom(atom_desc)
    try:
        return (
            DESC(ATOMIC_NUMBERS[element]),
            DESC(valence),
        )
    except KeyError:
        return (
            999,
            999,
        )

Cycle = NamedTuple('Cycle', [('i', int), ('n', int), ('j', int)])

def Small_Cycle(*args: Sequence[int]) -> Cycle:
    assert len(args) == 3, 'Cycles of length > 9 or containing neighbours with more than 9 bonds are not allowed (details: {0})'.format(args)
    return Cycle(*args)

GROUP_INDICES = (0, 1, 2, 3, 4)
LEFT_GROUP_INDEX, LEFT_ATOM_INDEX, RIGHT_ATOM_INDEX, RIGHT_GROUP_INDEX, CYCLES_INDEX = GROUP_INDICES

CHIRAL_MARKER = '*'

class Invalid_Dihedral_Angles(Exception):
    pass

class Dihedral_Fragment(object):
    def __init__(
        self,
        dihedral_string: Optional[str] = None,
        atom_list: Optional[Union[Tuple[List[str], str, str, List[str]], Tuple[List[str], str, str, List[str], str]]] = None,
        dihedral_angles: Optional[Tuple[List[float], List[float]]] = None,
    ) -> None:
        assert dihedral_string is not None or atom_list is not None, [dihedral_string, atom_list]

        if dihedral_string is not None:
            splitted_string = split_group_str(dihedral_string)
            neighbours_1 = [atom.upper() for atom in split_neighbour_str(splitted_string[LEFT_GROUP_INDEX])]
            self.atom_2 = splitted_string[LEFT_ATOM_INDEX].upper()
            self.atom_3 = splitted_string[RIGHT_ATOM_INDEX].upper()
            neighbours_4 = [atom.upper() for atom in split_neighbour_str(splitted_string[RIGHT_GROUP_INDEX])]
            self.cycles = (
                [
                    Small_Cycle(*list(map(int, cycle_str)))
                    for cycle_str in
                    split_neighbour_str(splitted_string[CYCLES_INDEX])
                ]
                if len(splitted_string) == 5
                else []
            )
            keep_stereoscopic_centers = True
        else:
            if len(atom_list) == 4:
                neighbours_1, self.atom_2, self.atom_3, neighbours_4 = atom_list
                self.cycles = []
            elif len(atom_list) == 5:
                neighbours_1, self.atom_2, self.atom_3, neighbours_4, self.cycles = atom_list
                self.cycles = [Small_Cycle(*c) for c in self.cycles]
            else:
                raise Exception('Wrong length of atom_list: {0}'.format(atom_list))
            keep_stereoscopic_centers = False

        self.neighbours_1, self.neighbours_4 = map(deque, (neighbours_1, neighbours_4))

        canonical_rep = self.__canonical_rep__(keep_stereoscopic_centers, dihedral_angles=dihedral_angles)

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

    def order_cycles(self: Any) -> None:
        self.cycles.sort(
            key=lambda cycle: (cycle.n, cycle.i, cycle.j),
        )

    def canonise_cycles(self: Any) -> None:
        if len(self.cycles) in (0, 1):
            pass
        elif len(self.cycles) == 2:
            cycle_0, cycle_1 = self.cycles
            should_order_left = (self.neighbours_1[cycle_0.i] == self.neighbours_1[cycle_1.i])
            should_order_right = (self.neighbours_4[cycle_0.j] == self.neighbours_4[cycle_1.j])
            self.cycles = [
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
                        (self.neighbours_1, self.neighbours_4),
                    )
                ],
            )

            if ENFORCE_CANONICAL_TRICYLIC:
                assert all(should_order.values()), 'Only symmetric (all similar or all different) environments are allowed for N-cycles where N >= 3 (fragment={0}, cycles={1}).'.format(
                    str(self),
                    self.cycles,
                )

            self.cycles = [
                Cycle(i, n, j)
                for (i, n, j) in
                zip(
                    *list(map(
                        sorted,
                        (Is, Ns, Js),
                    ))
                )
            ]

            def cycle_key(cycle: Cycle) -> Tuple[int, int, int, int]:
                return(
                    self.neighbours_1[cycle.i],
                    len([1 for other_cycle in self.cycles if cycle.i == other_cycle.i]),
                    self.neighbours_4[cycle.j],
                    len([1 for other_cycle in self.cycles if cycle.j == other_cycle.j]),
                )

    def __canonical_rep__(self, keep_stereoscopic_centers: bool, dihedral_angles: Optional[Tuple[List[float], List[float]]]) -> Any:
        other = copy(self)

        if not keep_stereoscopic_centers:
            other.sort_neighbours_renumber_cycles(dihedral_angles)
        other.canonise_cycles()
        other.order_cycles()
        other.flip_fragment_if_necessary()

        return other

    def flip_fragment_if_necessary(self: Any) -> None:
        # Compare the two central atoms and put the heavier one on the left
        if on_asc_atomic_number_then_asc_valence(self.atom_2) <  on_asc_atomic_number_then_asc_valence(self.atom_3):
            should_reverse = True
        elif on_asc_atomic_number_then_asc_valence(self.atom_2) == on_asc_atomic_number_then_asc_valence(self.atom_3):
            should_reverse = False

            if len(self.neighbours_1) < len(self.neighbours_4):
                should_reverse = True
            elif len(self.neighbours_1) == len(self.neighbours_4):
                # If identical central atoms, and same number of neighbours on both ends, try to resolve ambiguity one neighbour at a time
                for (neighbour_1, neighbour_4) in zip(self.neighbours_1, self.neighbours_4):
                    if on_asc_atomic_number_then_asc_valence(neighbour_1) < on_asc_atomic_number_then_asc_valence(neighbour_4):
                        should_reverse = True
                        break
                    elif on_asc_atomic_number_then_asc_valence(neighbour_1) == on_asc_atomic_number_then_asc_valence(neighbour_4):
                        pass
                    else:
                        break
            else:
                should_reverse = False
        else:
            should_reverse = False

        if should_reverse:
            self.reverse_dihedral()

    def sort_neighbours_renumber_cycles(self: Any, dihedral_angles: Optional[List[float]]) -> None:
        '''Order each neighbour list by alphabetical order, reordering the (possible) cycles to match those changes.'''
        if dihedral_angles is not None:
            left_dihedral_angles, right_dihedral_angles = dihedral_angles
            assert len(left_dihedral_angles) == len(self.neighbours_1) and len(right_dihedral_angles) == len(self.neighbours_4), [left_dihedral_angles, self.neighbours_1, right_dihedral_angles, self.neighbours_4]
            if not all(-180.0 <= angle <= 180.0 for angle in left_dihedral_angles + right_dihedral_angles):
                raise Invalid_Dihedral_Angles([left_dihedral_angles + right_dihedral_angles])
        else:
            left_dihedral_angles, right_dihedral_angles = [0.0 for _ in self.neighbours_1], [0.0 for _ in self.neighbours_4]

        def sorted_neighbours_permutation_dict(neighbours: List[str], angles: List[str]) -> Tuple[Deque[str], Dict[int, int]]:
            get_neighbour = lambda item: item[1][0]
            on_dihedral_angle_then_desc_atomic_number_and_valence = lambda item: (item[1][1], on_desc_atomic_number_then_desc_valence(get_neighbour(item)))

            sorted_neighbour_items = list(
                sorted(
                    enumerate(zip(neighbours, angles)),
                    key=on_dihedral_angle_then_desc_atomic_number_and_valence,
                    reverse=False,
                )
            )

            neighbour_items_deque = deque(sorted_neighbour_items)

            if DEBUG:
                print(
                    list(zip(neighbours, angles))
                )

            best_items = sorted(
                [
                    rotated_deque(neighbour_items_deque, n)
                    for n in range(len(sorted_neighbour_items))
                ],
                key=lambda _neighbours: tuple(
                    [
                        on_asc_atomic_number_then_asc_valence(neighbour)[0]
                        for (i, (neighbour, angle)) in _neighbours
                    ],
                ),
                reverse=True,
            )[0]

            permutation_dict = {
                i: j
                for (j, (i, _)) in enumerate(best_items)
            }
            return (
                deque(map(get_neighbour, best_items)),
                permutation_dict,
            )

        self.neighbours_1, permutation_1 = sorted_neighbours_permutation_dict(self.neighbours_1, left_dihedral_angles)
        self.neighbours_4, permutation_4 = sorted_neighbours_permutation_dict(self.neighbours_4, right_dihedral_angles)
        self.cycles = [
            Cycle(permutation_1[neighbour_id_1], cycle_length, permutation_4[neighbour_id_4])
            for (neighbour_id_1, cycle_length, neighbour_id_4) in self.cycles
        ]

    def reverse_dihedral(self) -> None:
        self.atom_2, self.atom_3 = self.atom_3, self.atom_2
        self.neighbours_1, self.neighbours_4 = self.neighbours_4, self.neighbours_1
        self.cycles = [x[::-1] for x in self.cycles]

    def is_left_chiral(self) -> bool:
        return len(set(self.neighbours_1)) >= 3

    def is_right_chiral(self) -> bool:
        return len(set(self.neighbours_4)) >= 3

    def is_chiral_fragment(self) -> bool:
        return (self.is_left_chiral() or self.is_right_chiral())
