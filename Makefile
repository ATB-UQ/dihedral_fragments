PYTHONPATH = PYTHONPATH="/home/$${USER}/ATB:/home/$${USER}/ATB-Dependencies"

test:
	$(PYTHONPATH) python fragment_dihedral.py
