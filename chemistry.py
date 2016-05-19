from fragment_dihedrals.fragment_dihedral import FragmentDihedral

CHEMICAL_GROUPS = (
    # Hydrocarbons
    ('alkane', 'J{3}|C|C|J{3}'),
    ('alkene', 'J{2}|C|C|J{2}'),
    ('alkyne', 'J|C|C|J'),

    # Groups containing halogens
    ('monofluoro', 'F,!X,!X|C|Z|J{2-3}'),
    ('difluoro', 'F,F,!X|C|Z|J{2-3}'),
    ('trifluoro', 'F{3}|C|Z|J{2-3}'),
    ('chloro', 'CL,!X,!X|C|Z|J{2-3}'),
    ('dichloro', 'CL,CL,!X|C|Z|J{2-3}%'),
    ('trichloro', 'CL{3}|C|Z|J{2-3}'),
    ('bromo', 'BR,!X,!X|C|Z|J{2-3}'),
    ('dibromo', 'BR,BR,!X|C|Z|J{2-3}'),
    ('tribromo', 'BR{3}|C|Z|J{2-3}'),
    ('iodo', 'I,!X,!X|C|Z|J{2-3}'),
    ('diiodo', 'I,I,!X|C|Z|J{2-3}'),
    ('triiodo', 'I{3}|C|Z|J{2-3}'),

    # Groups containing oxygen
    ('alcohol I', 'J,H,H|C|O|H'),
    ('alcohol II', 'C,C,H|C|O|H'),
    ('alcohol III', 'C{3}|C|O|H'),
    ('ketone', '%|C|C|O,C'),
    ('aldehyde', 'J{2-3}|C|C|O,H'),
    ('acyl halide', '%|C|C|O,X'),
    ('carboxylic acid', '%,O|C|O|H'),
    ('ester', 'J,O|C|O|C'),
    ('ether', 'J{3}|C|O|C'),

    # Groups containing nitrogen
    ('amide', 'O,J|C|N|J,J'),
    ('amine I', 'J{3}|C|N|H,H'),
    ('amine II', 'J{3}|C|N|C,H'),
    ('amine III', 'J{3}|C|N|C{2}'),
    ('ammonium ion', '%|J|N|J{3}'),
    ('ketimine I', 'H|N|C|C{2}'),
    ('ketimine II', 'C|N|C|C{2}'),
    ('aldimine I', 'H|N|C|C,H'),
    ('aldimine II', 'C|N|C|C,H'),
#    ('imide', ''),
    ('azide', 'N|N|N|C'),
    ('azo', 'C|N|N|C'),
    ('cyanate', 'C|O|C|N'),
    ('isocyanate', 'C|N|C|O'),
    ('nitrate', 'C|O|N|O{2}'),
    ('nitrile', 'N|C|C|Z{2-3}'),
    ('isonitrile', 'C|N|C|Z{2-3}'),
    ('nitro', '%|C|N|O{2}'),
    ('nitroso', '%|C|N|O'),
#    ('pyridyl', ''),

    # Groups containing sulfur
    ('thiol', 'J+|C|S|H'),
    ('thioether', 'J+|C|S|C'),
    ('disulfide', 'C|S|S|C'),
    ('sulfinyl', '%|C|S|C,O'),
    ('sulfonyl', '%|C|S|C,O,O'),
    ('sulfino', 'C,O|S|O|H'),
    ('sulfo', 'C,O,O|S|O|H'),
    ('thiocyanate', 'C|S|C|N'),
    ('isothiocyanate', 'C|N|C|S'),
    ('carbonothioyl', '%|J|C|S,J'),

    # Groups containing phosphorus
    ('phosphino', '%|C|P|C{2}'),
    ('phosphono', '%|C|P|O{3}'),
    ('phosphate', 'C|O|P|O{3}'),

    # Groups containing boron
    ('borono', '%|C|B|O{2}'),
)

CHEMICAL_GROUPS = map(
    lambda (moiety, pattern): (moiety, str(FragmentDihedral(dihedral_string=pattern))),
    CHEMICAL_GROUPS,
)

if __name__ == '__main__':
    print CHEMICAL_GROUPS
