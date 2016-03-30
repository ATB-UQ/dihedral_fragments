from itertools import product, combinations_with_replacement
from fragment_dihedral import FragmentDihedral

MONOVALENT = (1,)

ATOM_VALENCES = {
    'C': (2, 3, 4),
    'N': (2, 3, 4),
    'P': (3, 4),
    'O': (1, 2),
    'S': (1, 2),
    'H': MONOVALENT,
    'CL': MONOVALENT,
    'I': MONOVALENT,
    'BR': MONOVALENT,
    'F': MONOVALENT,
}

def is_monovalent(atom):
    return ATOM_VALENCES[atom] == MONOVALENT

def number_neighbours(atom):
    return [valence - 1 for valence in ATOM_VALENCES[atom]]

ATOMS = ATOM_VALENCES.keys()
CORE_ATOMS = [atom for atom in ATOMS if not is_monovalent(atom)]

FORBIDDEN_BONDS = (
    ('P', 'P'),
    ('C', 'P'),
    ('N', 'P'),
    ('BR', 'N'),
    ('CL', 'N'),
    ('I', 'N'),
    ('F', 'N'),
)

FORBIDDEN_BONDS = [sorted(bond) for bond in FORBIDDEN_BONDS]

def is_forbidden_bond(bond):
    return (sorted(bond) in FORBIDDEN_BONDS)

def main():
    for (atom_2, atom_3) in combinations_with_replacement(CORE_ATOMS, 2):

        if is_forbidden_bond((atom_2, atom_3)):
            continue

        for (len_neighbours_1, len_neighbours_4) in product(number_neighbours(atom_2), number_neighbours(atom_3)):
            neighbours_1 = combinations_with_replacement(ATOMS, len_neighbours_1)
            neighbours_4 = combinations_with_replacement(ATOMS, len_neighbours_4)
            for a, b in product(neighbours_1, neighbours_4):
                print FragmentDihedral(atom_list=(list(a), atom_2, atom_3, list(b)))

if __name__ == '__main__':
    main()
