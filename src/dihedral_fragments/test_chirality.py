from atb_api import API
from dihedral_fragments.capping import uncapped_molecule_for_dihedral_fragment
from atb_outputs.mol_data import MolData
from algorithm.atb.outputs import Output

if __name__ == '__main__':
    TEST_FRAGMENTS = ['S,O,N|P|O|C', 'S,N,O|P|O|C', 'O,H,C|C|C|O,H,C', 'H,H,H|C|C|H,H,H', 'CL,H|C|C|CL,H', 'C,C|N|C|O,H,C', 'C,C|N|C|O,C,H']

    for test_fragment in TEST_FRAGMENTS:
        print(test_fragment)
        uncapped_molecule = uncapped_molecule_for_dihedral_fragment(test_fragment)
        pdb = uncapped_molecule.dummy_pdb()

        mol_data = MolData(pdb)
        mol_data.completed = lambda x: False
        output = Output(mol_data, None, None)
        (_, dihedral_fragments), *_ = output.dihedral_fragments()
        try:
            assert dihedral_fragments == [test_fragment], (dihedral_fragments, [test_fragment])
        except AssertionError:
            print('ERROR: Missing fragment in uncapped molecule')
            print(pdb)
            raise

        uncapped_molecule.get_best_capped_molecule_with_ILP(enforce_octet_rule=True)
        pdb = uncapped_molecule.energy_minimised_pdb()

        mol_data = MolData(pdb)
        mol_data.completed = lambda x: False
        output = Output(mol_data, None, None)
        (_, dihedral_fragments), *_ = output.dihedral_fragments()
        try:
            assert test_fragment in dihedral_fragments, (test_fragment, dihedral_fragments)
        except AssertionError:
            print('ERROR: Missing fragment in capped molecule')
            print(pdb)
            raise

        print()
