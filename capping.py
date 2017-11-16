from typing import Any, List, Optional, Tuple
from pprint import pprint
from sys import stderr
from math import cos, sin, pi, sqrt

from fragment_capping.helpers.molecule import Uncapped_Molecule, Molecule
from fragment_capping.helpers.types_helpers import Fragment, Atom

from dihedral_fragments.dihedral_fragment import element_valence_for_atom, NO_VALENCE

def best_capped_molecule_for_dihedral_fragment(fragment_str: Fragment, debug: bool = False) -> Molecule:
    molecule = uncapped_molecule_for_dihedral_fragment(fragment_str).get_best_capped_molecule_with_ILP(
        debug=stderr if debug else None,
        enforce_octet_rule=True,
    )

    if debug:
        print(molecule)

    return molecule

def uncapped_molecule_for_dihedral_fragment(dihedral_fragment: Fragment, debug: bool = False) -> Uncapped_Molecule:
    if dihedral_fragment.count('|') == 3:
        neighbours_1, atom_2, atom_3, neighbours_4 = dihedral_fragment.split('|')
        cycles = []
        neighbours_1, neighbours_4 = neighbours_1.split(','), neighbours_4.split(',')
    elif dihedral_fragment.count('|') == 4:
        neighbours_1, atom_2, atom_3, neighbours_4, cycles = dihedral_fragment.split('|')
        neighbours_1, neighbours_4, cycles = neighbours_1.split(','), neighbours_4.split(','), cycles.split(',')
    else:
        raise Exception('Invalid dihedral_fragment: "{0}"'.format(dihedral_fragment))

    ids = [n for (n, _) in enumerate(neighbours_1 + [atom_2, atom_3] + neighbours_4, start=1)]

    neighbours_id_1, atom_id_2, atom_id_3, neighbours_id_4 = ids[0:len(neighbours_1)], ids[len(neighbours_1)], ids[len(neighbours_1) + 1], ids[len(neighbours_1) + 2:]
    CENTRAL_BOND = (atom_id_2, atom_id_3)

    elements = dict(
        list(zip(
            ids,
            [element_valence_for_atom(neighbour)[0] for neighbour in neighbours_1] + [atom_2, atom_3] + [element_valence_for_atom(neighbour)[0] for neighbour in neighbours_4],
        )),
    )

    valences = dict(
        list(zip(
            ids,
            [element_valence_for_atom(neighbour)[1] for neighbour in neighbours_1] + [len(neighbours_1) + 1, len(neighbours_4) + 1] + [element_valence_for_atom(neighbour)[1] for neighbour in neighbours_4],
        )),
    )

    bonds = (
        [
            (neighbour_id, atom_id_2)
            for neighbour_id in neighbours_id_1
        ]
        +
        [CENTRAL_BOND]
        +
        [
            (atom_id_3, neighbour_id)
            for neighbour_id in neighbours_id_4
        ]
    )

    def coordinates_for_atom_id(atom_id: int, d: float = 1.5) -> Tuple[float, float, float]:
        if atom_id == atom_id_2:
            return (-d / 2, 0, 0)
        elif atom_id == atom_id_3:
            return (d / 2, 0, 0)
        else:
            e = 0.5 * d
            f = sqrt(d ** 2 - e ** 2)
            assert e ** 2 + f ** 2 == d ** 2, (e ** 2 + f ** 2, d ** 2)
            if atom_id in neighbours_id_1:
                left_theta = 2 * pi / len(neighbours_id_1) * neighbours_id_1.index(atom_id)
                return (-d /2 - e, -f * cos(left_theta), f * sin(left_theta))
            elif atom_id in neighbours_id_4:
                right_theta = 2 * pi / len(neighbours_id_4) * neighbours_id_4.index(atom_id)
                return (d / 2 + e, f * cos(right_theta), f * sin(right_theta))
            else:
                raise Exception('Impossible id: {0}'.format(atom_id))

    molecule = Molecule(
        dict(
            list(zip(
                ids,
                [
                    Atom(
                        index=atom_id,
                        element=elements[atom_id],
                        valence=valences[atom_id],
                        capped=atom_id not in (neighbours_id_1 + neighbours_id_4),
                        coordinates=coordinates_for_atom_id(atom_id),
                    )
                    for atom_id in ids
                ],

            ))
        ),
        bonds,
        name=dihedral_fragment.replace('|', '_'),
    )

    if debug:
        print(molecule)

    for (i, n, j) in map(lambda cycle: map(int, cycle), cycles):
        i_id, j_id = neighbours_id_1[i], neighbours_id_4[j]
        if n == 0:
            # i and j are actually the same atoms
            del molecule.atoms[j_id]
            replace_j_by_i = lambda x: i_id if x == j_id else x
            molecule.bonds = {
                frozenset(map(replace_j_by_i, bond))
                for bond in molecule.bonds
            }
        else:
            NEW_ATOM_ID = -1
            NEW_ATOM = Atom(
                index=NEW_ATOM_ID, # This will get overwritten by Molecule.add_atom
                element='C',
                valence=NO_VALENCE,
                capped=False,
                coordinates=None,
            )
            atom_chain_id = [i_id] + [molecule.add_atom(NEW_ATOM) for i in range(n - 1)] + [j_id]
            new_bonds = zip(atom_chain_id[:-1], atom_chain_id[1:])
            molecule.add_bonds(new_bonds)

    if debug:
        print(molecule)

    return molecule


