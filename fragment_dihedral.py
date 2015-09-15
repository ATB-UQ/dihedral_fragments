from copy import deepcopy
class FragmentDihedral(object):

    def __init__(self, dihedral_string):
        splitted_string = dihedral_string.split("|")

        self.atom_1 = splitted_string[0][0]
        self.atom_2 = splitted_string[1][0]
        self.neighbours_2 = list(splitted_string[1].split('-')[1].split(','))
        self.atom_3 = splitted_string[2][0]
        self.neighbours_3 = list(splitted_string[2].split('-')[1].split(','))
        self.atom_4 = splitted_string[3][0]

    def __str__(self):
        return "{atom_1}|{atom_2}-{neighbours_2}|{atom_3}-{neighbours_3}|{atom_4}".format(
            atom_1=self.atom_1,
            atom_2=self.atom_2,
            neighbours_2=",".join(self.neighbours_2),
            atom_3=self.atom_3,
            neighbours_3=",".join(self.neighbours_3),
            atom_4=self.atom_4,
        )

    def __eq__(self, other):
        return self.__canonical_rep__().__dict__ == other.__canonical_rep__().__dict__

    def __ne__(self, other):
        return not self == other

    def __canonical_rep__(self):
        other = deepcopy(self)
        # Order each neighbour list by alphabetical order
        # WARNING: Keep in mind that maintaining order will be necessary for maintaining sterochemistry information
        other.neighbours_2.sort()
        other.neighbours_3.sort()
        # Compare the two neighbour list and put the first in lexical order on the leftmost side
        if [other.atom_3] + other.neighbours_3 <=  [other.atom_2] + other.neighbours_2:
            other.atom_1, other.atom_4 = other.atom_4, other.atom_1
            other.neighbours_2, other.neighbours3 = other.neighbours_3, other.neighbours_2
        return other

if __name__ == "__main__" :
    dihedral = FragmentDihedral("C|C-C,H|C-H,H|C")
    print dihedral
    dihedral_2 = FragmentDihedral("C|C-H,C|C-H,H|C")
    print dihedral_2.__canonical_rep__()
    print dihedral == dihedral_2
    dihedral_3 = FragmentDihedral("C|C-H,C|C-H,C|C")
    print dihedral_3 == dihedral_2

