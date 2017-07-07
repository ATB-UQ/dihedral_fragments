from typing import List, Any, Optional

from dihedral_fragments.dihedral_fragment import GROUP_SEPARATOR, on_asc_atomic_number_then_asc_valence, NEIGHBOUR_SEPARATOR

class Improper(object):
    def __init__(self, central: Optional[str] = None, neighbours: Optional[List[str]] = None, improper_str: Optional[str] = None) -> None:
        if improper_str is not None:
            fields = improper_str.split(GROUP_SEPARATOR)
            assert len(fields) == 2, fields
            self.central, self.neighbours = fields[0], fields[1].split(NEIGHBOUR_SEPARATOR)
        else:
            self.central, self.neighbours = central, neighbours

        self.neighbours = sorted(
            self.neighbours,
            key=on_asc_atomic_number_then_asc_valence,
            reverse=True,
        )

    def __str__(self) -> str:
        return 'Improper(central={0}, neighbours={1})'.format(self.central, self.neighbours)

if __name__ == '__main__':
    print(Improper(improper_str='C|C,C,H'))
