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

```
>>> from dihedral_fragments.dihedral_fragment import Dihedral_Fragment
>>> left_substituents, left_atom, right_atom, right_substituents = (['C', 'N', 'H'], 'N', 'C', ['H', 'H', 'H'])
>>> str(Dihedral_Fragment(atom_list=(left_substituents, left_atom, right_atom, right_substituents), dihedral_angles=([0., 120., -120.], [-120, 0., 120])))
'N,H,C|N|C|H,H,H'
>>> str(Dihedral_Fragment(atom_list=(left_substituents, left_atom, right_atom, right_substituents), dihedral_angles=([0., -120., 120.], [-120, 0., 120])))
'N,C,H|N|C|H,H,H'
>>> str(Dihedral_Fragment(atom_list=(left_substituents, left_atom, right_atom, right_substituents)))
'N,C,H|N|C|H,H,H'
```

### Canonise a fragment

```
>>> from dihedral_fragments.dihedral_fragment import Dihedral_Fragment
>>> fragment = 'C,N,O|C|N|H,H,H'
>>> canonical_fragment = str(Dihedral_Fragment(fragment))
>>> canonical_fragment
'H,H,H|N|C|O,C,N'
```
