# VASP Test
Instructions
1. add `module load vasp/5.4.4` to your `.bashrc` file and reload it or reconnect 
2. run file `00_make_molecule.py` by running `python3 00_make_molecule.py`
3. verify that the file `start.traj` was created 
4. run the second file with slurm by typing `sbatch 01_submit_job.py` (you might need a __init__.py file for kultools.py to be recognized as a module). 
5. Check the `job.out` and `job.err` files for errors
