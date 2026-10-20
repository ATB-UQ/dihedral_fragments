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
cd /home/atb/ATB/dihedral_fragments/src/dihedral_fragments && python -m pytest -q   # 10 pass, 2 xfail (known bugs)
```

`test.py` is mostly print-only; `test_grammar.py` holds the real assertions (canonical-form
invariances, round trips, pattern matching) and two strict xfails recording known bugs (review
2026-10: two-letter elements with a valence digit such as `CL1` mis-sort, and
`sql_pattern_matching_for` interpolates user text into SQL). The `Makefile` targets still point at
the retired Python 3.5 / `ATB_ONE` layout.

## What still works, and what does not

| module | state |
|---|---|
| `dihedral_fragment.py` (`Dihedral_Fragment`, `canonical_representation_for`), `improper.py`, `regex.py`, `atomic_numbers.py`, `pattern_matching.py`, `chemistry.py` | live: used by `core` (`core/atb/_outputs/dihedral_fragments.py`), `atb_bonded_motifs.keys`, `dihedral_gap_service.keys`, and the website (`pattern_matching`, `chemistry`: restricted tab and API call `molids_with_dihedral_fragment`) |
| `capping.py`, `tag_predictor.py`, `fragment_generator.py` | import, but nothing imports them (hand-run only; `tag_predictor` compiles patterns and prints at import) |

Removed in the 2026-10 review (they raised ImportError on `atb_api.HTTPError`, used a dead 2017
host, and had no importer): `cap.py`, `molecule_for_fragment.py`, `tasks/cap_all.py`,
`latex_table.py`, `test_capping.py`, `test_chirality.py`; see git history.

## Dependencies

`setup.cfg` declares only `cairosvg` (now unused); `pattern_matching` needs `jinja2` (undeclared),
`capping` needs `fragment_capping`, `tag_predictor` needs `atb_client`; all come from the
shared venv. `dihedral_fragment.py` itself needs only the standard library.

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
`website`; uses `fragment_capping` (capping), `atb_outputs` and `atb_client` (the broken
modules above). No environment variables, no vendor binaries, no DB.

# Citation / Attribution

To cite this work, please use the following [Zenodo DOI](https://zenodo.org/badge/latestdoi/95523757).
