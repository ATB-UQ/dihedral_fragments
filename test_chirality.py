from atb_api import API
from dihedral_fragments.capping import uncapped_molecule_for_dihedral_fragment
from atb_outputs.mol_data import MolData
from algorithm.atb.outputs import Output

if __name__ == '__main__':
    TEST_FRAGMENTS = ['O,H,C|C|C|O,H,C', 'H,H,H|C|C|H,H,H', 'H,CL|C|C|CL,H']

    for test_fragment in TEST_FRAGMENTS:
        uncapped_molecule = uncapped_molecule_for_dihedral_fragment(test_fragment)
        pdb = uncapped_molecule.dummy_pdb()

        mol_data = MolData(pdb)
        mol_data.completed = lambda x: False
        output = Output(mol_data, None, None)
        (_, dihedral_fragments), *_ = output.dihedral_fragments()
        try:
            assert dihedral_fragments == [test_fragment], (dihedral_fragments, [test_fragment])
        except AssertionError:
            print(pdb)
            raise
