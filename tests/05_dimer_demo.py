#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#SBATCH -p high
#SBATCH --job-name=opt_try
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64 -c 2
#SiBATCH --exclusive

"""
Created on Thu Mar 26 17:10:56 2020
@author: ark245
"""
import os, sys
sys.path.insert(0,'/home/ark245/lib/onboarding_DFT/')
from ase import io
from kul_tools import KulTools as KT

kt = KT(gamma_only=False,structure_type='zeo',calculation_type='dimer')
atoms = io.read('dimer_start.traj')
kt.set_structure(atoms)
kt.set_overall_vasp_params({'gga':'RP'})
kt.run()
