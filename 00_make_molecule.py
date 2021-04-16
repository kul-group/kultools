#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np  
from ase.build import molecule
from ase.visualize import view 
from ase import view


atoms = molecule('H2O')
atoms.set_cell(8*np.identity(3))
atoms.center()
view(atoms)
io.write('water.traj')



