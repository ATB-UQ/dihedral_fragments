from argparse import ArgumentParser
from typing import Any, List, Optional

from dihedral_fragments.molecule_for_fragment import molid_after_capping_fragment, Fragment, Too_Many_Permutations, Molecule_Not_In_ATB, PDB_Structure_Not_Found, ATB_Molecule_Running

def parse_args():
    parser = ArgumentParser()

    parser.add_argument('--fragment', type=Fragment, help='')
    parser.add_argument('--profile', action='store_true', help='')
    parser.add_argument('--debug', action='store_true', help='Run overly-verbose debugging')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()
    if not args.profile:
        print(
            molid_after_capping_fragment(args.fragment, quick_run=False, debug=args.debug),
        )
    else:
        from cProfile import runctx as profile_run
        from pstats import Stats

        time_file='.time'
        profile_run(
            'molid_after_capping_fragment(fragment, quick_run=True)',
            {},
            dict(molid_after_capping_fragment=molid_after_capping_fragment, fragment=args.fragment, debug=args.debug),
            filename=time_file,
        )

        stats = Stats(time_file).sort_stats('tottime', 'cumtime')
        stats.print_stats(500)
