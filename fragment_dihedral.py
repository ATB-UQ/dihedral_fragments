from copy import deepcopy
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

    def __str__(self):
        return "{neighbours_1}|{atom_2}|{atom_3}|{neighbours_4}".format(
            neighbours_1=','.join(self.neighbours_1),
            atom_2=self.atom_2,
            atom_3=self.atom_3,
            neighbours_4=','.join(self.neighbours_4),
        )

    def __eq__(self, other):
        return self.__canonical_rep__().__dict__ == other.__canonical_rep__().__dict__

    def __ne__(self, other):
        return not self == other

    def __canonical_rep__(self):
        other = deepcopy(self)
        # Order each neighbour list by alphabetical order
        # WARNING: Keep in mind that maintaining order will be necessary for maintaining sterochemistry information
        other.neighbours_1.sort()
        other.neighbours_4.sort()
        # Compare the two neighbour list and put the first in lexical order on the leftmost side
        if [other.atom_2] + other.neighbours_1 >=  [other.atom_3] + other.neighbours_4:
            other.atom_2, other.atom_3 = other.atom_3, other.atom_2
            other.neighbours_1, other.neighbours_4 = other.neighbours_4, other.neighbours_1
        return other

if __name__ == "__main__" :
    dihedral_1 = FragmentDihedral("C,C,H|C|C|C,H,H")
    print dihedral_1
    dihedral_2 = FragmentDihedral("C,H,H|C|C|C,C,H")
    print dihedral_2
    print dihedral_1 == dihedral_2
    dihedral_3 = FragmentDihedral(atom_list=(['H', 'H'], 'C', 'C', ['Cl', 'Cl']))
    print dihedral_3
    print dihedral_3 == dihedral_2

