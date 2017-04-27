from dihedral_fragments.dihedral_fragment import Dihedral_Fragment, sql_pattern_matching_for, re_pattern_matching_for

def test_canonical_rep() -> None:
    test_cases = (
        ('H,C4,H|SI|C|C2,H,C4', 'C4,H,H|SI|C|C4,C2,H'), # M(SI) > M(C)
        ('H,C4,H|C|C|C4,H', 'C4,H,H|C|C|C4,H'), #M(C) == M(C), 3 > 2
        ('H,C4,H|C|C|C3,H,H', 'C4,H,H|C|C|C3,H,H'), #M(C) == M(C), 3 ==3, 4 > 3
        ('O2|C|C|N2,H', 'N2,H|C|C|O2'), # M(O) > M(N)
        ('F,C4|C|C|O1,N3', 'F,C4|C|C|O1,N3'), # M(F) > M(O)
    )

    for (dihedral_string, solution) in test_cases:
        answer = str(Dihedral_Fragment(dihedral_string))
        assert answer == solution, '"{0}" (answer) != "{1}" (expected)'.format(answer, solution)

        cyclic_answer = str(Dihedral_Fragment(answer))
        assert cyclic_answer == answer, '"{0}" (answer) != "{1}" (expected)'.format(cyclic_answer, answer)

def test_patterns() -> None:
    test_cases = (
        'X,X|C|C|X,X',
        '!X+|N|C|CX',
        'C,X,B,D|N|N|H,_,C',
        'C+|CC|C|C+',
        'CL,CL,X{2}|C|Z|CX',
        'J+|C|S|H',
        '!J{2-4}|C|C|C',
        'X{2},I|C|Z|%',
        '!X{2},I|C|Z|%',
        'CL{2-3}|C|C|BR{3-5}',
        'J,J|C|C|H',
        'H|C|C|J{3}',
    )

    for pattern in test_cases:
        print(pattern)
        print(sql_pattern_matching_for(pattern))
        print()

def test_cyclic_fragments() -> None:
    cyclic_fragment = Dihedral_Fragment(atom_list=(['H', 'H', 'C'], 'C', 'C', ['H', 'C', 'H'], [[2, 0, 1]]))
    print(str(cyclic_fragment))
    assert str(cyclic_fragment) == 'C,H,H|C|C|C,H,H|000', cyclic_fragment

    polycyclic_fragment, answer = Dihedral_Fragment(atom_list=(['H', 'C', 'C'], 'C', 'C', ['C', 'C', 'H'], [[2, 2, 1], [1, 3, 0]])), 'C,C,H|C|C|C,C,H|020,131'
    print(str(polycyclic_fragment))
    assert str(polycyclic_fragment) == answer, '{0} != {1} (expected)'.format(
        str(polycyclic_fragment),
        answer,
    )
    print(Dihedral_Fragment(str(polycyclic_fragment)))

    polycyclic_fragment, answer = Dihedral_Fragment(atom_list=(['H', 'N', 'O'], 'C', 'C', ['C', 'C', 'H'], [[2, 2, 1], [1, 3, 0]])), 'C,C,H|C|C|O,N,H|020,131'
    print(str(polycyclic_fragment))
    assert str(polycyclic_fragment) == answer, '{0} != {1} (expected)'.format(
        str(polycyclic_fragment),
        answer,
    )
    print(Dihedral_Fragment(str(polycyclic_fragment)))

    assert str(Dihedral_Fragment('C,C,C|C|C|C,C,C|002,101,200')) == 'C,C,C|C|C|C,C,C|000,101,202'

    try:
        Dihedral_Fragment('C,C,N|C|C|C,C,C|002,101,200')
        raise Exception('This should have failed.')
    except AssertionError:
        print('Non-symmetric N=3 rings failed as expected.')

    try:
        Dihedral_Fragment('C,C,C|C|C|C,C,C|0100')
        raise Exception('This should have failed.')
    except AssertionError:
        print('Rings length > 9 failed as expected.')

def test_atom_list_init() -> None:
    fragment = Dihedral_Fragment(atom_list=(['C'], 'C', 'C', ['C']))

def test_misc() -> None:
    dihedral_1 = Dihedral_Fragment("C,C,H|C|C|C,H,H")
    print(dihedral_1)
    dihedral_2 = Dihedral_Fragment("C,H,H|C|C|C,C,H")
    print(dihedral_2)
    print(dihedral_1 == dihedral_2)
    dihedral_3 = Dihedral_Fragment(atom_list=(['H', 'H'], 'C', 'C', ['Cl', 'Cl']))
    print(dihedral_3)
    print(dihedral_3 == dihedral_2)

def test_chiral_str() -> None:
    dihedral_1 = Dihedral_Fragment("C,C,H|C|C|C,H,H")
    print(dihedral_1.__str__())
    print(dihedral_1.__str__(flag_chiral_sides=True))

if __name__ == "__main__" :
    #test_atom_list_init()
    test_patterns()
    test_canonical_rep()
    test_chiral_str()
    test_cyclic_fragments()
    test_misc()

    assert re_pattern_matching_for('Z,%|Z|Z|Z,%', debug=True)('C,H|C|C|C,H') == True
    assert re_pattern_matching_for('Z|Z|Z|Z,%', debug=True)('C,H|C|C|C,H') == False
    assert re_pattern_matching_for('Z|Z|Z|Z,%', debug=True)('C|C|C|C,H') == True
    assert re_pattern_matching_for('N,J|C|C|J{3}', debug=True)('N,H|C|C|C,H,H') == True
    assert re_pattern_matching_for('N,J|C|C|J{2}', debug=True)('N,H|C|C|C,H') == True

    print(re_pattern_matching_for('Z|Z|Z|Z', debug=True))

    print(sql_pattern_matching_for('J{3}|C|C|J{3}'))
