from copy import deepcopy
class FragmentDihedral(object):

    def __init__(self, dihedral_string):
        splitted_string = dihedral_string.split("|")

        self.atom1 = splitted_string[0][0]
        self.neighbours1 = list(splitted_string[0][1:])
        self.atom2 = splitted_string[1][0]
        self.atom3 = splitted_string[1][1]
        self.atom4 = splitted_string[2][0]
        self.neighbours4 = list(splitted_string[2][1:])

    def __str__(self):
        return "{atom1}{neighbours1}|{atom2}{atom3}|{atom4}{neighbours4}".format(atom1=self.atom1, atom2=self.atom2, atom3=self.atom3, atom4=self.atom4, neighbours1="".join(self.neighbours1), neighbours4="".join(self.neighbours4))

    def __eq__(self, other):
        return self.__canonical_rep__().__dict__ == other.__canonical_rep__().__dict__

    def __ne__(self, other):
        return not self == other

    def __canonical_rep__(self):
        other = deepcopy(self)
        # Order each neighbour list by alphabetical order
        # WARNING: Keep in mind that maintaining order will be necessary for maintaining sterochemistry information
        other.neighbours1.sort()
        other.neighbours4.sort()
        # Compare the two neighbour list and put the first in lexical order on the leftmost side
        if [ other.atom4 ] + other.neighbours4 <= [ other.atom1 ] + other.neighbours1:
            other.atom1, other.atom4 = other.atom4, other.atom1
            other.neighbours1, other.neighbours4 = other.neighbours4, other.neighbours1
        return other

if __name__ == "__main__" :
    dihedral = FragmentDihedral("CAB|CC|CHH")
    print dihedral
    dihedral2 = FragmentDihedral("CHH|CC|CAB")
    print dihedral2.__canonical_rep__()
    print dihedral == dihedral2
