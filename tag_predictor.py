from sys import stderr

from API_client.api import API
from dihedral_fragments.fragment_dihedral import re_pattern_matching_for
from dihedral_fragments.chemistry import CHEMICAL_GROUPS

DEBUG = True

CHEMICAL_GROUPS_MATCHING_PATTERNS = [
    (moiety, re_pattern_matching_for(pattern, debug=DEBUG, metadata=moiety))
    for (moiety, pattern) in CHEMICAL_GROUPS if pattern
]

def dihedrals(molecule):
    return molecule.dihedral_fragments

def tags_for_dihedral(dihedral_string):
    tags = [moiety for (moiety, matching_function) in CHEMICAL_GROUPS_MATCHING_PATTERNS if matching_function(dihedral_string)]
    assert len(tags) <= 1, 'No dihedral ({0}) should be matched by more than one rule: {1}'.format(dihedral_string, tags)
    return tags

def tags_for_molecule(molecule):
    return set(
        reduce(
            lambda acc, e: acc + e,
            [tags_for_dihedral(dihedral) for dihedral in dihedrals(molecule)],
            [],
        ),
    )

def yes_or_no(query):
    choice = eval(input(query))
    if choice == 'y':
       return True
    elif choice == 'n':
       return False
    else:
       stderr.write("Please respond with 'yes' or 'no'")

IGNORE_FILE = '.ignore'

def get_ignored_molids():
    from os.path import exists
    if exists(IGNORE_FILE):
        with open(IGNORE_FILE) as fh:
            return set([int(line) for line in fh.read().splitlines()])
    else:
        return set()

if __name__ == '__main__':
    assert tags_for_dihedral('CL,C,H|C|C|H,H,H') == ['chloro']
    #assert tags_for_dihedral('C|N|C|C,H') == ['']
    assert tags_for_dihedral('CL,CL,H|C|C|H,H,H') == ['dichloro']
    assert tags_for_dihedral('C,C|N|C|H,H,H') == ['amine III']
    assert tags_for_dihedral('H,H,H|C|C|O,C') == ['ketone']
    assert tags_for_dihedral("O,H|C|C|C,H") == ['aldehyde']
    assert tags_for_dihedral('H,H,H|C|C|CL,CL,CL') == ['trichloro']
    assert tags_for_dihedral('O,O,O|P|O|C') == ['phosphate']

    ignored_molids = get_ignored_molids()

    fh = open(IGNORE_FILE, 'a')

    DATA_SETS_TAGS = set(['Shivakumar et al.', 'Marenich et al.', 'Mobley et al.', 'SAMPL0', 'SAMPL1', 'SAMPL2', 'SAMPL4'])

    api = API(debug=True, api_format='pickle')
    molecules = api.Molecules.search(tag='Mobley et al.', max_atoms=15)
    for molecule in molecules:
        if molecule.molid in ignored_molids:
            continue

        print(molecule)
        new_tags = tags_for_molecule(molecule)
        old_tags = set(molecule.tags) - DATA_SETS_TAGS

        if old_tags != new_tags:
            print('Old tags are: {0}'.format(list(old_tags)))
            print('New tags are: {0}'.format(list(new_tags)))

            try:
                if yes_or_no('Add those ({0}) tags ?'.format(new_tags - old_tags)):
                    print('Adding new tags ...')
                    for tag_name in new_tags - old_tags:
                        print(molecule.tag(tag_name=tag_name))
                    print('Successfully added new tags')
                else:
                    fh.write(str(molecule.molid) + '\n')
            except KeyboardInterrupt:
                fh.close()
                raise
        print()
