PYTHONPATH = PYTHONPATH="/home/$${USER}/ATB_ONE"

test:
	$(PYTHONPATH) python3 test.py

fragment_generator:
	$(PYTHONPATH) python3 fragment_generator.py

tables:
	$(PYTHONPATH) python3 latex_table.py

errors:
	/usr/local/python35/bin/pylint -E $$(find . -name '*.py')
.PHONY: errors

mypy:
	MYPYPATH=/home/$$USER/ATB_ONE /usr/local/python35/bin/mypy dihedral_fragment.py
.PHOMY:
