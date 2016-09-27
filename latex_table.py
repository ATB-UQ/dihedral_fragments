from collections import Iterable
from jinja2 import Template

from .chemistry import CHEMICAL_GROUPS

def dihedral(pattern):
    if pattern:
        return '\\dihedral{' + protect(pattern) + '}'
    else:
        return ''

def protect(pattern):
    return pattern.replace('%', '\\%').replace('{', '\\{').replace('}', '\\}')

def indent_str(latex_str, number_indents):
    return '\n'.join(['\t'*number_indents + line for line in latex_str.splitlines()])

def latex_render_table(rows, row_formatters=None, indent=4, number_columns=1):
    assert len(rows) >= 1

    row_length = len(rows[0])

    if not row_formatters:
        row_formatters = [lambda x:str(x)]*row_length

    if len(rows) % number_columns == 0:
        extra_row = None
    else:
        padding_columns = len(rows) % number_columns
        try:
            empty_element = ['' for field in rows[0]]
        except TypeError:
            empty_element = ['']

        extra_row = reduce(
            lambda acc, e: acc + e,
            [tuple(empty_element)]*(number_columns - padding_columns) + list(rows[-padding_columns:]),
            (),
        )

    boundary = len(rows) // number_columns

    extended_rows = [[] for i in range(boundary)]

    for (i, row) in enumerate(rows):
        extended_rows[i % boundary] += row

    table = Template('''
\\begin{tabular}{{ tabular_align }}
{%- for row in rows %}
\t{{ process_row(row) }} \\\\
{%- endfor %}
\\end{tabular}
''').render(
        rows=extended_rows + ([extra_row] if extra_row else []),
        dihedral=dihedral,
        process_row=lambda row: ' & '.join([formatter(field) for (field, formatter) in zip(row, row_formatters*number_columns)]),
        tabular_align='{' + '|'.join(['c' for x in range(row_length)]*number_columns) + '}'
    )
    return indent_str(table, indent)

def latex_caption(caption, indent=4):
    return indent_str(r'\caption{' + caption + '}', indent)

if __name__ == "__main__" :
    print(latex_render_table(
        rows=[(name, pattern) for (name, pattern) in CHEMICAL_GROUPS if pattern],
        row_formatters=(lambda x: x, lambda x: dihedral(x)),
        number_columns=3,
    ))

    from urllib.request import urlopen
    import json
    response = urlopen('https://atb.uq.edu.au/index.py?filter_0=40&tab=dihedral_fragments&format=json&items_per_page=30')
    top30 = json.loads(response.read())['main']

    print(''.join([
        latex_caption(''),
        latex_render_table(
            rows=top30,
            row_formatters=(lambda x: dihedral(str(x)), lambda x: str(x)),
            number_columns=3,
        ),
    ]))
