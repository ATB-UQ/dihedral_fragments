from typing import List

REGEX_START_ANCHOR, REGEX_END_ANCHOR, REGEX_OR_OPERATOR = ('^', '$', '|')

BACKSLASH = '\\'
COMMA = ','
ESCAPED_COMMA = ';'
UNESCAPE_COMMA = lambda x: x.replace(ESCAPED_COMMA, COMMA)
ATOM_CHARACTERS, MAX_ATOM_CHARACTERS = 'A-Z', 2
VALENCE_CHARACTERS, MAX_VALENCE_CHARACTERS = '0-9', 1
ONE_ATOM = '[' + ATOM_CHARACTERS + VALENCE_CHARACTERS + ']' + '{1,' + str(MAX_ATOM_CHARACTERS + MAX_VALENCE_CHARACTERS) + '}'
ONE_NUMBER = '[0-9]'
ANY_NUMBER_OF_ATOMS = '[' + ATOM_CHARACTERS + VALENCE_CHARACTERS + ESCAPED_COMMA + ']*'
ONE_OR_MORE_TIMES = '+'
GROUP = lambda n: r'\{n}'.format(n=n)

def REGEX_NOT_SET(pattern: str) -> str:
    return '[^(' + str(pattern) + ')]'

def REGEX_GROUP(index: int) -> str:
    return BACKSLASH + str(index)

def REGEX_SET(*args: List[str]) -> str:
    return '[' + ''.join([x for x in args]) + ']'

def REGEX_OR(*args: List[str]) -> str:
    return '(' + REGEX_OR_OPERATOR.join(args) + ')'

def REGEX_ESCAPE(pattern: str, flavour: str = 'sql') -> str:
    if flavour == 'sql':
        return BACKSLASH + str(pattern)
        return BACKSLASH + BACKSLASH + str(pattern)
    elif flavour == 're':
        return '[' + str(pattern) + ']'
    else:
        raise Exception()

def REGEX_AT_LEAST(pattern: str, escape_plus: bool = True) -> str:
    return str(pattern) + (BACKSLASH if escape_plus else '') + '+'

def FORMAT_ESCAPED(pattern: str) -> str:
    return pattern.replace('{', '{{').replace('}', '}}')

def FORMAT_UNESCAPED(pattern: str) -> str:
    return pattern.replace('{{', '{').replace('}}', '}')

def NOT(pattern: str) -> str:
    return '!' + str(pattern)

def CAPTURE(pattern: str) -> str:
    return '(' + str(pattern) + ')'

def ESCAPE(pattern: str) -> str:
    return BACKSLASH + str(pattern)

def exactly_N_times_operator(what: str = ONE_ATOM, how_many_times: int = ONE_NUMBER) -> str:
    return CAPTURE(what) + ESCAPE('{') + CAPTURE(how_many_times) + ESCAPE('}')

def N_to_M_times_operator(what: str = ONE_ATOM, N: int = ONE_NUMBER, M: int = ONE_NUMBER) -> str:
    return CAPTURE(ONE_ATOM) + ESCAPE('{') + CAPTURE(ONE_NUMBER) + ESCAPE('-') + CAPTURE(ONE_NUMBER) + ESCAPE('}')


