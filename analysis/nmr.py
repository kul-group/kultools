#!/usr/bin/env python
from ase.io import read, write
from glob import glob
import sys
from ase.neighborlist import NeighborList, natural_cutoffs

kwd = " (absolute, valence and core) "  # id's the line where relevant data is.

# User specifies what folders to loop through
path_list = sys.argv[1:]
path_list.sort()
print(path_list)

for path in path_list:
    # Read in atoms object, find index of N next to Sn
    # atoms object indexed starting at 0, OUTCAR index starts at 1; thus the +1
    atoms = read(path + "/../input_for_nmr.traj")
    nl = NeighborList(natural_cutoffs(atoms), bothways=True, self_interaction=False)
    nl.update(atoms)

    sn_index = [i.index for i in atoms if i.symbol == "Sn"][0]
    neigh_ind = nl.get_neighbors(sn_index)[0]
    ind = [i for i in neigh_ind if atoms[i].symbol == "N"][0] + 1

    # Read the OUTCAR and report ISO_SHIFT for the atom
    with open(path + "/OUTCAR", "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if kwd in line:
                data_line_ind = i + ind
                break  # takes only the first matching section

        iso_shift = lines[data_line_ind].split()[
            4
        ]  # Gets the 5th column containing iso_shift (incl. G=0 contribution)
        print(iso_shift)
