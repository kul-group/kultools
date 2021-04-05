#!/usr/bin/env python
#SBATCH -J md -q regular -N 4 -t 47:00:00 -C knl --output=job.out --error=job.err --ntasks-per-node=64

import shutil
from ase.parallel import *
from ase import io
from ase.atoms import np
import os,sys
import time
import shutil
from ase.optimize import BFGS,FIRE, QuasiNewton
from glob import glob
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.optimize.gpmin import gpmin
path_home = os.environ['HOME']
sys.path.insert(0,path_home+'/lib/zeolite_tools')
from identify_zeolite_atoms import analyze_zeolite_atoms

path_home = os.environ['HOME']
if path_home.startswith('/global/homes'):
	host_name = 'cori'
elif path_home.startswith('/home/'):
	host_name = 'hpc1'
elif path_home.startswith('/Users/'):
	host_name = 'local'
elif path_home.startswith('/home1/'):
	host_name = 'stampede'

if 'hpc1' in host_name:
    os.environ['VASP_PP_PATH']='/home/ark245/programs/vasp5.4.4/pseudopotentials/pseudo54/'
    os.environ['VASP_COMMAND']='module load vasp/5.4.4pl2-vtst; NTASKS=`echo $SLURM_TASKS_PER_NODE|tr \'(\' \' \'|awk \'{print $1}\'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo " $NTASKS * $NNODES " | bc`; echo "num_cpu=" $NCPU; srun -n $NCPU vasp_gam| tee -a op.vasp'
elif 'cori' in host_name or 'edison' in host_name or 'nid' in host_name:
	print(os.environ['HOSTNAME'])
	os.environ['VASP_PP_PATH']='/global/homes/a/ark245/pseudopotentials/PBE54/'
	os.environ['VASP_COMMAND']='NTASKS=`echo $SLURM_TASKS_PER_NODE|tr \'(\' \' \'|awk \'{print $1}\'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo " $NTASKS * $NNODES " | bc`; echo "num_cpu=" $NCPU; srun -n $NCPU vasp_gam | tee -a op.vasp'

def analyze_zeolite_structure(atoms):
	total_T_atoms = sum([1 for a in atoms if a.symbol in ['Si','Al','Zr','Hf','Sn']])
	total_O_zeolite_atoms = total_T_atoms*2
	total_zeolite_atoms = total_T_atoms + total_O_zeolite_atoms
	print("Total zeolite atoms =", total_zeolite_atoms)
	return total_zeolite_atoms

def initialize_magmoms(atoms, total_zeolite_atoms):
	magmoms=atoms.get_initial_magnetic_moments()
	for a in atoms:
		if a.symbol in ['Cu', 'Ni', 'Co', 'Fe']:
			magmoms[a.index] = 2.5
		elif a.symbol in ('Al','Zr','Hf','Sn'):
			magmoms[a.index] = 0.5
		elif a.symbol == 'N':
			magmoms[a.index] = 0.25
		elif a.index > total_zeolite_atoms and a.symbol == 'O':
			magmoms[a.index] = 2.0
		else:
			magmoms[a.index] = 0.0
	atoms.set_initial_magnetic_moments(magmoms)
	#print (magmoms)
	return atoms

def assign_calculator(atoms,my_encut=400):
	calc_for_opt = Vasp(kpts=(1,1,1),
		potim=1.0,
		encut=my_encut,
		ispin=2,
		nsw=50,
		prec='Normal',
		istart=1,
		isif=2,
		ismear=0,
		sigma=0.05,
		nelmin=4,
		nelmdl=-4,
		nwrite=1,
		icharg=2,
		lasph=True,
		ediff=1E-6,
		ediffg=-0.05,
		ibrion=2,
		lcharg=False,
		lwave=False,
		laechg=False,
		voskown=1,
		algo='Fast',
		lplane=True,
		lreal='Auto',
		isym=0,
		xc='PBE',
		lorbit=10,
		nupdown=-1,
		npar=4,
		nsim=4,
		ivdw=12)
	atoms.set_calculator(calc_for_opt)
	return (atoms)

def opt(atoms,dir_name='opt',**kwargs):

	total_zeolite_atoms = analyze_zeolite_structure(atoms)
	atoms = assign_calculator(atoms) # Using ASE calculator
	atoms.calc.set(ibrion=2, potim = 0.5,**kwargs)
	#for key,val in kwargs.items():
	#	print(key,val)
	#	atoms.calc.set(key=val) #lreal=False) #kwargs[key])
	atoms = initialize_magmoms(atoms,total_zeolite_atoms)


	cwd = os.getcwd()
	encut_for_dir = atoms.calc.float_params['encut']
	directory = dir_name + '_' + str(encut_for_dir)
	if not os.path.exists(directory):
		os.makedirs(directory)
	os.chdir(directory)

	energy = atoms.get_potential_energy() # Run vasp here

	#change back
	os.chdir(cwd)

	return (atoms)

def dimer(atoms,**kwargs):
	total_zeolite_atoms = analyze_zeolite_structure(atoms)
	atoms = assign_calculator(atoms) # Using ASE calculator
	atoms = initialize_magmoms(atoms,total_zeolite_atoms)
	atoms.calc.set(ibrion=3, potim=0, ichain=2, iopt=7, ddr = 0.01, drotmax = 2, **kwargs) #, nupdown = my_nupdown)
	atoms.get_potential_energy()
	print ('Done!  \nGetting energies.')

def opt_fine(atoms,**kwargs):
	total_zeolite_atoms = analyze_zeolite_structure(atoms)
	atoms = assign_calculator(atoms) # Using ASE calculator
	atoms = initialize_magmoms(atoms,total_zeolite_atoms)
	atoms.calc.set(ibrion=1, potim = 0.05, nsw = 5, ediffg=-0.03, **kwargs)
	energy = atoms.get_potential_energy() # Run vasp here
	return (atoms)

def vib(atoms,**kwargs):
	total_zeolite_atoms = analyze_zeolite_structure(atoms)
	atoms = assign_calculator(atoms) # Using ASE calculator
	atoms = initialize_magmoms(atoms,total_zeolite_atoms)
	atoms.calc.set(ibrion=5,potim=0.02,nsw=1,lwave=False,lcharg=False, **kwargs)
	energy = atoms.get_potential_energy()
	return(atoms)

def md(atoms,dir_name='md',**kwargs):
	total_zeolite_atoms = analyze_zeolite_structure(atoms)
	atoms = assign_calculator(atoms,my_encut=300) # Using ASE calculator
	atoms = initialize_magmoms(atoms,total_zeolite_atoms)
	atoms.calc.set(ibrion=0,potim=0.5,tebeg=298, nsw=100000, isif=2, smass=0, **kwargs)

	cwd = os.getcwd()
	encut_for_dir = atoms.calc.float_params['encut']
	directory = dir_name + '_' + str(encut_for_dir)
	if not os.path.exists(directory):
		os.makedirs(directory)
	os.chdir(directory)

	energy = atoms.get_potential_energy() # Run vasp here

	#change back
	os.chdir(cwd)
	return(atoms)

def apply_constraints(atoms,constraints_type):
	#atoms = io.read('vasprun.xml') #temp_opt.traj')
	if constraints_type == 'zeo':
		indices, count, indices_atomtype, count_atomtype, atomtype_list = analyze_zeolite_atoms(atoms)
		indices_to_fix = []
		for i, atom_type in enumerate(atomtype_list):
			if 'framework-Si' in atom_type or \
				 'framework-O' in atom_type or  \
					'framework-Al' in atom_type :
				indices_to_fix.append(i)
		c = FixAtoms(indices=[a.index for a in atoms if a.index in indices_to_fix])
		atoms.set_constraint(c)
	else:
		pass
	return (atoms)

def vib_workflow(atoms,type_vib='std', **kwargs):
	if 'ts' not in type_vib and 'isolated' not in type_vib:
		assert len(atoms.constraints) == 1
		cwd = os.getcwd()
		directory='opt_fine'
		if not os.path.exists(directory):
			os.makedirs(directory)
		os.chdir(directory)
		atoms = opt_fine(atoms, **kwargs)
		os.chdir(cwd)

	if 'ts' in type_vib:
		directory = 'vib_ts'
		assert len(atoms.constraints) == 1
	elif 'isolated' in type_vib:
		directory='vib_isolated'
	else:
		directory = 'vib'
	if not os.path.exists(directory):
		os.makedirs(directory)

	cwd = os.getcwd()
	os.chdir(directory)
	atoms = vib(atoms, **kwargs)
	os.chdird(cwd)
	return (atoms)

def run_mode(mode,atoms,constraints_type,**kwargs):
	# initial opt
	if mode == 'opt':
		atoms = opt(atoms,**kwargs)

	if 'solv' in mode:
		if not 'solv-opt' in mode and not 'solv-spe' in mode:
			print('ERROR: Check mode')
			sys.exit()
		cwd = os.getcwd()
		if 'opt' in mode:
			dir_name = 'opt'
			my_nsw = 100
		elif 'spe' in mode:
			dir_name = 'spe'
			my_nsw = 0
		atoms = opt(atoms,dir_name=dir_name,lwave=True,lsol=False,nsw=my_nsw)

		directory = mode #solv-spe or solv-opt
		if not os.path.exists(directory):
			os.makedirs(directory)
		os.chdir(directory)
		shutil.copyfile('../spe/WAVECAR','WAVECAR')
		atoms = opt(atoms,dir_name='solv-spe',nsw=0,lwave=False,lsol=True)
		os.chdir(cwd)

	if mode == 'md':
		md(atoms)

	# vibrations
	if 'vib' in mode  and 'ts' not in mode and 'isolated' not in mode:
		atoms = apply_constraints(atoms,constraints_type)
		atoms = vib_workflow(atoms,type_vib='std',**kwargs)
	elif 'vib' in mode and 'ts' in mode:
		atoms = apply_constraints(atoms,constraints_type)
		atoms = vib_workflow(atoms,type_vib='ts',**kwargs)
	elif 'vib' in mode and 'isolated' in mode:
		atoms = vib_workflow(atoms,type_vib='isolated', **kwargs)

	# dimer
	if mode == 'dimer':
		if not os.path.exists('NEWMODECAR'):
			atoms = io.read('dimer_start.traj')
			#atoms = apply_constraints(atoms)
		else:
			atoms = io.read('CENTCAR',format='vasp')
			shutil.copyfile('NEWMODECAR', 'MODECAR')
			shutil.copyfile('OUTCAR', 'OUTCAR.bak')
			shutil.copyfile('vasprun.xml', 'vasprun.xml.bak')

			#atoms = apply_constraints(atoms)
		dimer(atoms,**kwargs)

	return atoms

##

atoms = io.read('manual_start.traj')
mode = 'opt'
calc_args = {'encut':300,'xc':'PBE','kpts':(1,1,1),'kpar':1,'nsw':20}
atoms = run_mode(mode,atoms,constraints_type='none',**calc_args)

mode = 'md'
print('MD starting')
calc_args = {'encut':300,'xc':'PBE','kpts':(1,1,1),'kpar':1,'nsw':1000000}
atoms = run_mode(mode,atoms,constraints_type='none',**calc_args)
print('MD done')
