#!/usr/bin/env python

##SBATCH --exclude=agate-0,agate-1
#SBATCH --nodes=1 --partition=high
#SBATCH --ntasks-per-node=64
#SBATCH --ntasks-per-core=1
#SBATCH --threads-per-core=1
#SBATCH --output=job.out
#SBATCH --error=job.err
#SBATCH --time=28:00:00
#SBATCH --verbose

import shutil
from ase.parallel import *
from ase.io import read, write
from ase import io
from ase.atoms import np
import os
import time
import shutil
from ase.optimize import BFGS,FIRE, QuasiNewton
from glob import glob
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.optimize.gpmin import gpmin

os.environ['VASP_PP_PATH']='/home/sours/programs/vasp_PP'
os.environ["MPIRUN_OPTIONS"] = "--verbose --map-by core --display-map --report-bindings "
os.environ["MPIRUN_OPTIONS"] += "--mca btl_openib_allow_ib true --mca btl_openib_if_include mlx5_0:1 "
os.environ["MPIRUN_OPTIONS"] += "--mca btl_openib_warn_nonexistent_if 0 "
os.environ["MPIRUN_OPTIONS"] += "--mca btl_openib_warn_no_device_params_found 0 "
os.environ["MPIRUN_OPTIONS"] += "--mca pml ob1 --mca btl openib,self,vader "
os.environ["MPIRUN_OPTIONS"] += "--mca mpi_cuda_support 0 "
os.environ['VASP_COMMAND'] = "mpirun $MPIRUN_OPTIONS -np $SLURM_NTASKS vasp_gam"
