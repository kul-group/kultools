#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 17:10:56 2020

@author: ark245
"""
import os
os.chdir('/Users/ark245/Box/codebase/dev_kul_tools/')
from ase import io
from kul_tools import KulTools as KT

kt = KT(gamma_only=False,structure_type='zeo')
kt.set_calculation_type('opt')
atoms = io.read('manual_start.traj')
#traj = io.read('vasprun.xml',':')
kt.set_structure(atoms)
kt.set_overall_vasp_params({'gga':'RP','nupdown':5})
kt.run()

#
print(kt.calc_default.int_params['ismear'])
print(kt.calc_default.special_params['lreal'])
print(kt.calculation_type)


# Calculation type 
# opt --> needs atoms 
# opt_fine --> needs constraints 
# vib --> needs constraints 
# dimer --> needs path annd atoms 
# neb --> needs traj 

