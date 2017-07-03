# Dihedral Fragments v1.0

Author: Bertrand Caron

## Requirements

The root `dihedral_fragments` directory must be in your `PYTHONPATH` shell variable.
Failure to do so will result in an `ImportError` being raised:

```
$ PYTHONPATH="" python3 -c 'from dihedral_fragments.dihedral_fragment import Dihedral_Fragment'
Traceback (most recent call last):
  File "<string>", line 1, in <module>
ImportError: No module named 'dihedral_fragments'
```

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
