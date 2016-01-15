from jinja2 import Template
from fragment_dihedral import CHEMICAL_GROUPS

def dihedral(pattern):
    return '\\dihedral{' + protect(pattern) + '}'

def protect(pattern):
    return pattern.replace('%', '\\%')

def latex_render_table(rows, row_formatters=None, indent=4):
    assert len(rows) >= 1

    row_length = len(rows[0])

    if not row_formatters:
        row_formatters = [lambda x:str(x)]*row_length

    table = Template('''
\\begin{tabular}{{ tabular_align }}
{%- for row in rows %}
\t{{ process_row(row) }} \\\\
{%- endfor %}
\\end{tabular}
''').render(
        rows=rows,
        dihedral=dihedral,
        process_row=lambda row: ' & '.join([formatter(field) for (field, formatter) in zip(row, row_formatters)]),
        tabular_align='{' + '|'.join(['c' for x in range(row_length)]) + '}'
    )
    indented_table = '\n'.join(['\t'*indent + line for line in table.splitlines()])
    return indented_table

if __name__ == "__main__" :
    print latex_render_table(
        rows=CHEMICAL_GROUPS,
        row_formatters=(lambda x: x, lambda x: dihedral(x)),
    )
