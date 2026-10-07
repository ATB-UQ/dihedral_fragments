[![DOI](https://zenodo.org/badge/95523757.svg)](https://zenodo.org/badge/latestdoi/95523757)

# Dihedral Fragments v1.0

Author: Bertrand Caron

A **dihedral fragment** is the string key for a rotatable bond and its first neighbours,
e.g. `C,H,H|C|C|H,H,H` (left substituents | left atom | right atom | right substituents),
canonicalised so that equivalent environments compare equal and, optionally, carrying the
stereochemistry of the substituent order. This package defines the class and its grammar, and a
few tools built on it (regex/SQL pattern matching, capping a fragment into a real molecule).

**Status: live-support.** It sits at the head of the dihedral subsystem: *fragment keys* (here) ->
`dihedral_gap_service` (finds fragments with no archetype) -> `dihedral_scan_service` (QM scans)
-> archetype fit and publication (`atb_graph_helpers`) -> the IFP/archetype library `core`
uses. Only the core of the package is alive; see "What still works" below.

## Install and test

A normal editable install in the shared venv (no `PYTHONPATH` needed; the old note about putting
the repo root on `PYTHONPATH` is obsolete):

```bash
uv pip install --python /home/atb/ATB/.venv/bin/python -e /home/atb/ATB/dihedral_fragments --no-deps
cd /home/atb/ATB/dihedral_fragments/src/dihedral_fragments && python -m pytest -q test.py   # 7 pass
```

`test_chirality.py` and `test_capping.py` currently fail at collection. The `Makefile` targets
(`make test`, `mypy`) still point at the retired Python 3.5 / `ATB_ONE` layout.

## What still works, and what does not

| module | state |
|---|---|
| `dihedral_fragment.py` (`Dihedral_Fragment`, `canonical_representation_for`), `improper.py`, `regex.py`, `atomic_numbers.py`, `pattern_matching.py` | live: used by `core` (`core/atb/_outputs/dihedral_fragments.py`, `_outputs/mtb.py`, `helpers/fragment_dihedrals.py`), `atb_bonded_motifs.keys`, `dihedral_gap_service.keys`, and the website's `molecules` API (`pattern_matching`) |
| `capping.py` | imports; wraps `fragment_capping` to cap a fragment |
| `cap.py`, `molecule_for_fragment.py`, `tag_predictor.py`, `fragment_generator.py`, `tasks/cap_all.py` | **broken in the current venv**: they do `from atb_api import API, HTTPError, ATB_Mol` (the modern `atb_api_public` no longer exports `HTTPError`) and use a dead 2017 host (`scmb-atb.biosci.uq.edu.au/atb-uqbcaron`). Used only by two hand-run website tasks (`cap_all_ATB_fragments.py`, `submit_protein_fragments.py`) |

## Dependencies

`setup.cfg` declares only `cairosvg`; `fragment_capping`, `atb_outputs` and `atb_api_public` are
commented out of it (they were `git+ssh` URLs) and come from the shared venv. Importing the
live modules needs only the standard library.

## Use cases

### Create new fragment

Fragment is defined by 5 variables: `left_substituents`, `left_atom`, `right_atom`, `right_substituents` and `cycles`.
For the example fragment
`(['A1','A2','A3'], 'B', 'C', ['D1', 'D2'], [])`,
`dihedral_angles` should be a tuple of two lists of angles:

```([d(A1, B, C, DX), d(A2, B, C, DX), d(A3, B, C, DX)], [d(D1, C, B, AX), d(D2, C, B, AX)])```

where `AX` denotes one atom in `{'A1', 'A2', 'A3'}`, `DX` one atom in `{'D1', 'D2'}` and `d()` is the dihedral angle function (in either radians or degrees, units do not alter the relative ordering).
If no `dihedral_angles` are provided, the fragment is canonised (substituents ordered by descending atomic number, then descending number of bonded partners).
If they are provided, the stereochemistry is encoded into the fragment (order of the susbstituents).

```
>>> from dihedral_fragments.dihedral_fragment import Dihedral_Fragment
>>> left_substituents, left_atom, right_atom, right_substituents, cycles = (['C', 'N', 'H'], 'N', 'C', ['H', 'H', 'H'], [])
>>> str(Dihedral_Fragment(atom_list=(left_substituents, left_atom, right_atom, right_substituents, cycles)))
'N,C,H|N|C|H,H,H'
>>> str(Dihedral_Fragment(atom_list=(left_substituents, left_atom, right_atom, right_substituents, cycles), dihedral_angles=([0., 120., -120.], [-120, 0., 120])))
'N,H,C|N|C|H,H,H'
>>> str(Dihedral_Fragment(atom_list=(left_substituents, left_atom, right_atom, right_substituents, cycles), dihedral_angles=([0., -120., 120.], [-120, 0., 120])))
'N,C,H|N|C|H,H,H'
```

### Canonise a fragment

Note how the stereochemistry of the fragment is conserved.

```
>>> from dihedral_fragments.dihedral_fragment import Dihedral_Fragment
>>> fragment = 'C,N,O|C|N|H,H,H'
>>> str(Dihedral_Fragment(fragment))
'H,H,H|N|C|C,N,O'
```

## Platform wiring

Used by `core`, `atb_bonded_motifs`, `dihedral_gap_service`, `dihedral_parameterization` and
`website`; uses `fragment_capping` (capping), `atb_outputs` and `atb_api_public` (the broken
modules above). No environment variables, no vendor binaries, no DB.

# Citation / Attribution

To cite this work, please use the following [Zenodo DOI](https://zenodo.org/badge/latestdoi/95523757).
