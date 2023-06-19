from contextlib import redirect_stdout

from dihedral_fragments.molecule_for_fragment import molid_after_capping_fragment

missing_fragments = {'O,C,H|C|C|O,H,H', 'C|O|C|C,H,H', 'H|S|C|C,H,H', 'P|O|C|C,H,H', 'H,H|N|C|N,C', 'O,C,H|C|C|N,C,H', 'H,H|N|C|O,N', 'H|O|C|C,C,H', 'O,N,H|C|C|C,H,H', 'C,H,H|C|C|C,H', 'H|O|C|C,C', 'H|O|C|O,C', 'O,H,H|C|C|N,C,H', 'C,C|N|C|O,C,H', 'C,H,H|C|C|C,H,H', 'H|N|C|N,N', 'H,H,H|N|C|C,H,H', 'S,H,H|C|C|N,C,H', 'H,H|N|C|N,N', 'C,C,H|C|C|C,H,H', 'C,H|N|C|C,H,H', 'O,C,H|C|C|C,C,H', 'H,H|N|C|O,C', 'O,O,O|P|O|C', 'C|N|C|C,H,H', 'H,H|N|C|C,H,H', 'N,H,H|C|C|C,H,H', 'H|O|C|C,H,H', 'S,H,H|C|C|C,H,H', 'C,H|N|C|N,N', 'C|O|C|C,C,H', 'O,C,H|C|C|O,C,H', 'C|S|C|C,H,H', 'C,H|C|C|C,H', 'C,H,H|C|C|O,N', 'H|O|C|C', 'C,H,H|C|C|N,C', 'C|O|C|N,C,H', 'N,C,H|C|C|C,C,H', 'O,C,H|C|C|C,H,H', 'O,H,H|C|C|N,H,H', 'N,C,H|C|C|C,H,H', 'C,H,H|C|C|O,O', 'O,C,H|C|C|N,H,H', 'C,C,C|N|C|C,H,H', 'C|O|C|O,C', 'O,N,H|C|C|O,C,H', 'C,H,H|C|C|C,C'}

if __name__ == '__main__':
    for fragment in missing_fragments:
        try:
            with redirect_stdout(None):
                molids = molid_after_capping_fragment(fragment, quick_run=False)
            print(fragment, molids)
        except:
            print(fragment)
            raise
