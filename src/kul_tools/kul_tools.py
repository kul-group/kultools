# Kulkarn group @ UC Davis 
"""Definition of the Zeotype class.

This module defines the central object in analyzing Zeotypes 
"""
import sys
import os
import shutil
import signal

from datetime import datetime
import numpy as np 

from ase import Atoms
from ase.neighborlist import NeighborList
try:
    from ase.neighborlist import natural_cutoffs
except:
    from ase.utils import natural_cutoffs
from ase.visualize import view
from ase.calculators.vasp import Vasp

from pymatgen.io.vasp.outputs import Vasprun

class KulTools:
    """KulTools class that provides all the necessary tools for running simulations. Currently targetted towards using vasp. """
    def __init__(self,gamma_only=False,structure_type=None,calculation_type='spe',structure=None, is_stop_eligible=False): 
        """
        """

        """_summary_

        Args:
            gamma_only (bool, optional): _description_. Defaults to False.
            structure_type (str, optional): one of 'zeo','mof','metal','gas-phase','insulators', 'semiconductors'. Defaults to None.
            calculation_type (str, optional): _description_. Defaults to 'spe'.
            structure (_type_, optional): _description_. Defaults to None.
            is_stop_eligible(bool, False): 
        """

        self.working_dir = os.getcwd()
        
        self.hpc = self.identify_hpc_cluster()
        print('KT: HPC= %s' %self.hpc)
        self.gamma_only = gamma_only
        print('KT: VASP_GAMMA= %s' %self.gamma_only)
        
        self.structure_type = structure_type
        self.calculation_type = calculation_type
        self.main_dir = os.getcwd()
        self.structure = structure 
        self.is_stop_eligible = is_stop_eligible

        self.identify_vasp_eviron()
        print('KT: VASP_PP_PATH= %s' %self.vasp_pp_path)
        print('KT: VASP_COMMAND= %s' %self.vasp_command)

        self.assign_default_calculator()

    def identify_hpc_cluster(self):
        path_home = os.environ['HOME']
        if path_home.startswith('/global/homes'):
            if os.path.exists('/global/project/projectdirs'):
                host_name = 'cori'
            else:
                host_name = 'perlmutter'
        elif path_home.startswith('/home/'):
            if os.environ['SLURM_SUBMIT_HOST']=='hpc1':
                host_name = 'hpc1'
            elif os.environ['SLURM_SUBMIT_HOST']=='hpc2':
                host_name = 'hpc2'
        elif path_home.startswith('/Users/'):
            host_name = 'local'
        elif path_home.startswith('/home1/'):
            host_name = 'stampede'
        elif path_home.startswith("/jet/home/"):
            host_name = "bridges2"
        elif path_home.startswith("/g/g91/"):
            host_name = "quartz"
        else:
            print('Check cluster settings')
            sys.exit()
        return host_name
            
    def identify_vasp_eviron(self):
        if self.hpc == 'hpc1':
            os.environ['VASP_PP_PATH']='/home/ark245/programs/vasp5.4.4/pseudopotentials/pseudo54'
            if self.gamma_only: 
                vasp_exe = 'vasp_gam'
            else:
                vasp_exe = 'vasp_std'
            os.environ['VASP_COMMAND']='module load vasp/5.4.4pl2-vtst; NTASKS=`echo $SLURM_TASKS_PER_NODE|tr \'(\' \' \'|awk \'{print $1}\'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo " $NTASKS * $NNODES " | bc`; echo "num_cpu=" $NCPU; srun -n $NCPU %s | tee -a op.vasp' % vasp_exe

        elif self.hpc == 'hpc2':
            os.environ['VASP_PP_PATH']= '/home/sours/programs/vasp_PP' # '/home/ark245/programs/pseudopotentials/pseudo54' 
            if self.gamma_only:
                vasp_exe = 'vasp_gam'
            else:
                vasp_exe = 'vasp_std'
            os.environ['VASP_COMMAND']='NTASKS=`echo $SLURM_TASKS_PER_NODE|tr \'(\' \' \'|awk \'{print $1}\'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo " $NTASKS * $NNODES " | bc`; echo "num_cpu=" $NCPU; $(which mpirun) --map-by core --display-map --report-bindings --mca btl_openib_allow_ib true --mca btl_openib_if_include mlx5_0:1 --mca btl_openib_warn_nonexistent_if 0 --mca btl_openib_warn_no_device_params_found 0 --mca pml ob1 --mca btl openib,self,vader --mca mpi_cuda_support 0 -np $NCPU %s | tee -a op.vasp' % vasp_exe 
            print(os.environ['VASP_COMMAND'])

        elif self.hpc == 'cori':
            os.environ['VASP_PP_PATH']='/global/homes/a/ark245/pseudopotentials/PBE54'
            if self.gamma_only: 
                vasp_exe = 'vasp_gam'
            else:
                vasp_exe = 'vasp_std'
            os.environ['VASP_COMMAND']='NTASKS=`echo $SLURM_TASKS_PER_NODE|tr \'(\' \' \'|awk \'{print $1}\'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo " $NTASKS * $NNODES " | bc`; echo "num_cpu=" $NCPU; srun -n $NCPU %s | tee -a op.vasp' % vasp_exe
        elif self.hpc == 'stampede':
            os.environ['VASP_PP_PATH']='/home1/05364/ark245/pseudopotentials/PBE54'
            if self.gamma_only: 
                vasp_exe = 'vasp_gam_vtst'
            else:
                vasp_exe = 'vasp_std_vtst'
            os.environ['VASP_COMMAND']='module load vasp/5.4.4; export OMP_NUM_THREADS=1;rm op.vasp; mpirun -np $SLURM_NTASKS %s | tee op.vasp' % vasp_exe
        elif self.hpc == "bridges2":
            print(os.environ["HOSTNAME"])
            os.environ["VASP_PP_PATH"] = "/jet/home/rgoel/uo2/vasp_PP"
            os.environ[
                "ASE_VASP_COMMAND"
            ] = 'mpirun -np $SLURM_NTASKS /opt/packages/VASP/VASP5/PGI/vasp_std'
        elif self.hpc == 'local':
            os.environ['VASP_PP_PATH']='local_vasp_pp'
            if self.gamma_only: 
                vasp_exe = 'vasp_gam'
            else:
                vasp_exe = 'vasp_std'
            os.environ['VASP_COMMAND']='local_%s' % vasp_exe
        else: 
            if (os.environ.get('VASP_PP_PATH', None) is None) or ((os.environ.get('VASP_COMMAND', None) or os.environ.get('ASE_VASP_COMMAND', None)) is None):
                print('Check cluster settings')
                sys.exit()
        self.vasp_pp_path = os.environ['VASP_PP_PATH']
        self.vasp_command = os.environ.get('VASP_COMMAND', None) or os.environ.get('ASE_VASP_COMMAND', None)

    def assign_default_calculator(self):
        """Sets a default calculator regadless of the structure type"""
        self.calc_default = Vasp(kpts=(1,1,1),
            potim=0.5,
            encut=500,
            ispin=2,
            nsw=50,
            prec="Normal",
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
            lorbit=11,
            nupdown=-1,
            npar=4,
            kpar=1,
            nsim=4,
            ivdw=12)
        self.modify_calc_according_to_structure_type()
        
    def set_overall_vasp_params(self,overall_vasp_params):
        self.overall_vasp_params = overall_vasp_params
        
    def modify_calc_according_to_structure_type(self):
        if self.structure_type == 'zeo' or self.structure_type == 'mof': 
            pass
        elif self.structure_type == 'metal': 
            self.calc_default.set(sigma=0.2,ismear=1)
        elif self.structure_type in ['insulators', 'semiconductors']: 
            # ref: https://www.vasp.at/wiki/index.php/ISMEAR#Summary
            self.calc_default.set(sigma=0.05,ismear=-5)
        elif self.structure_type == 'gas-phase': 
            self.calc_default.set(kpts=(1,1,1),lreal=False)
        else:
            raise ValueError("Unknown structure_type = %s" % self.structure_type)
            
    def set_calculation_type(self,calculation_type):
        assert self.calculation_type in ['spe','opt','opt_fine'], "Unknown calculation_type = %s" % self.calculation_type
        self.calculation_type = calculation_type
        
    def set_structure(self,atoms_or_traj):
        if isinstance(atoms_or_traj, Atoms):
            atoms = atoms_or_traj
            self.structure = atoms 
            self.allowed_calculation_types = ['opt','opt_fine','vib']
            self.structure_istraj = False
        if isinstance(atoms_or_traj, list) and isinstance(atoms_or_traj[0], Atoms): 
            #print('Dealing with an traj object')
            traj = atoms_or_traj
            self.structure = traj
            self.allowed_calculation_types = ['neb']
            self.structure_istraj = True
            
    def get_structure(self):
        return self.structure
        
    def _change_to_dir(self,dir_name):
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        os.chdir(dir_name)

    def run_dft(self,atoms,dir_name):
        if self.is_stop_eligible:
            signal.signal(signal.SIGUSR1, self.checkpoint)
        atoms.set_calculator(self.calc)
        atoms.calc.set(**self.overall_vasp_params)
        #if self.calculation_type == 'opt' or self.calculation_type == 'vib':
        encut_for_dir = atoms.calc.float_params['encut']
        kpts_for_dir = ''.join([str(val) for val in atoms.calc.input_params['kpts']])
        func_for_dir = atoms.calc.input_params['xc']
        directory = dir_name + '_' + str(func_for_dir) + '_' + str(encut_for_dir) + '_' + str(kpts_for_dir)

        self._change_to_dir(directory)
        energy = atoms.get_potential_energy() # Run vasp here
        os.chdir(self.working_dir) 

        # Check for convergence after optimization
        #f = atoms.get_forces()
        #forces=np.sqrt(np.square(f[:,0]) +np.square(f[:,1]) + np.square(f[:,2]))
        #max_forces = max(forces)
        #magmoms = atoms.get_magnetic_moments()
        #mag_oszi = atoms.get_magnetic_moment()
        return (atoms)
    
    def run(self):
        if self.calculation_type == 'opt':
            self.structure_after = self.run_opt()#atoms,dir_name)
        elif self.calculation_type == 'opt_fine':
            self.structure_after = self.run_opt_fine()
        elif self.calculation_type == 'vib':
            self.structure_after = self.run_vib()
        elif self.calculation_type == 'solv':
            self.structure_after = self.run_solv()
        elif self.calculation_type == 'md':
            self.structure_after = self.run_md()
        elif self.calculation_type == 'spe':
            self.structure_after = self.run_spe()
        return self.structure_after

    def run_spe(self):
        """Runs a simple single point energy"""
        atoms = self.structure
        dir_name = 'spe' 
        self.calc = self.calc_default
        self.calc.set(nsw=0)
        new_atoms = self.run_dft(atoms,dir_name)
        return new_atoms
    
    def run_opt(self):
        """Runs a simple optimization"""
        atoms = self.structure
        self.calc = self.calc_default
        dir_name = 'opt'
        new_atoms = self.run_dft(atoms,dir_name)
        return new_atoms
    
    def run_opt_fine(self):
        """Runs a finer optimization"""
        atoms = self.structure
        dir_name = 'opt_fine' 
        self.calc = self.calc_default
        self.calc.set(ibrion=1, potim = 0.05, nsw = 50, ediffg=-0.03)
        new_atoms = self.run_dft(atoms,dir_name)
        return new_atoms
    
    def run_vib(self):
        """Runs a simple vib calculation"""
        atoms = self.structure
        dir_name = 'vib' 
        self.calc = self.calc_default
        self.calc.set(ibrion=5,potim=0.02,nsw=1)
        new_atoms = self.run_dft(atoms,dir_name)
        return new_atoms
    
    def run_solv(self,lrho=False):
        """Runs a simple solvation calculation"""
        atoms = self.structure
        dir_name = 'spe' 
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        self.calc = self.calc_default
        self.calc.set(potim=0.0,nsw=5,lwave=True,lsol=False,prec='Accurate')
        new_atoms = self.run_dft(atoms,dir_name)
        
        atoms = new_atoms
        dir_name = 'solv-spe'
        if not os.path.exists(dir_name):
            os.makedirs(dir_name)
        self.calc = self.calc_default
        shutil.copyfile('spe/WAVECAR','solv-spe/WAVECAR')
        self.calc.set(potim=0.0,nsw=3,lwave=True,lsol=True,prec='Accurate')
        new_atoms = self.run_dft(atoms,dir_name)
        
        if lrho: 
            atoms = new_atoms
            dir_name = 'solv-spe-rho'
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
            self.calc = self.calc_default
            shutil.copyfile('spe-spe/WAVECAR','solv-spe-rho/WAVECAR')
            self.calc.set(potim=0.0,nsw=0,lwave=True,lsol=True,prec='Accurate',lrhob=True,lrhoion=True)
            new_atoms = self.run_dft(atoms,dir_name)
        
        return new_atoms
        
    def run_md(self):
        """Runs a finer optimization"""
        atoms = self.structure
        dir_name = 'md'
        self.calc = self.calc_default
        self.calc.set(nsw=100000,ibrion=0,tebeg=298, isif=2, smass=0)
        new_atoms = self.run_dft(atoms,dir_name)
        return new_atoms
    
    
    def checkpoint(self, signum, _):
        print(f'Handling signal {signum} ({signal.Signals(signum).name}).')
        with open("STOPCAR", "w") as stopcar:
            stopcar.write("LSTOP = .TRUE.\n")
        print("Wrote stopcar and waiting for the run to end @",os.getcwd())

        # below code ref: https://stackoverflow.com/a/74112271/7630458
        # scontrol requeue doc:https://slurm.schedmd.com/scontrol.html#SECTION_COMMANDS
        # starttime slurm ref: https://slurm.schedmd.com/scontrol.html#SECTION_JOBS---SPECIFICATIONS-FOR-UPDATE-COMMAND
        # update command ref: https://stackoverflow.com/a/73239001/7630458
        # os.system('scontrol requeue $SLURM_JOB_ID')
        # os.system('scontrol update jobid=$SLURM_JOB_ID StartTime=now+1h')
        # os.system('sbatch --begin=now+120minutes -D={0} {1}'.format(root,sys.argv[0]))
        # os.system('sbatch --dependency=afterany:$SLURM_JOB_ID -D={0} {1}'.format(root,sys.argv[0]))



#    if not 'solv-opt' in mode and not 'solv-spe' in mode: 
#            print('ERROR: Check mode')
#            sys.exit()
#        cwd = os.getcwd()
#        if 'opt' in mode: 
#            dir_name = 'opt'
#            my_nsw = 100
#        elif 'spe' in mode: 
#            dir_name = 'spe'
#            my_nsw = 0
#        atoms = opt(atoms,dir_name=dir_name,lwave=True,lsol=False,nsw=my_nsw)
#
#        directory = mode #solv-spe or solv-opt
#        if not os.path.exists(directory):
#            os.makedirs(directory)
#        os.chdir(directory)
#        shutil.copyfile('../spe/WAVECAR','WAVECAR')
#        atoms = opt(atoms,dir_name='solv-spe',nsw=0,lwave=False,lsol=True)

    def run_specific_calcualtion_type():
        if self.calculation_type == 'opt':
            atoms = run_opt()
        elif self.calculation_type == 'opt_fine':
            atoms = run_opt_fine()
        elif self.calculation_type == 'vib':
            atoms = run_vib()
        elif self.calculation_type == 'solv':
            atoms = run_vib()
    
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
                
            
    
        # vibrations 
        if 'vib' in mode  and 'ts' not in mode and 'isolated' not in mode:
            atoms = apply_constraints(atoms)
            atoms = vib_workflow(atoms, **kwargs)
        elif 'vib' in mode and 'ts' in mode: 
            atoms = apply_constraints(atoms) 
            atoms = vib_workflow(atoms,type_vib='ts',**kwargs)
        elif 'vib' in mode and 'isolated' in mode: 
            atoms = vib_workflow(atoms,type_vib='isolated', **kwargs)
    
        # dimer 
        if mode == 'dimer':
            if not os.path.exists('NEWMODECAR'):
                atoms = io.read('dimer_start.traj')
                atoms = apply_constraints(atoms)
                dimer(atoms) 
            else: 
                atoms = io.read('CENTCAR',format='vasp')
                shutil.copyfile('NEWMODECAR', 'MODECAR') 
                shutil.copyfile('OUTCAR', 'OUTCAR.bak') 
                shutil.copyfile('vasprun.xml', 'vasprun.xml.bak') 
                
                atoms = apply_constraints(atoms)
                dimer(atoms)
        
        return atoms
        










    
