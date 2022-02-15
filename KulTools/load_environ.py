import sys, os


def identify_hpc_cluster():
    path_home = os.environ["HOME"]
    if path_home.startswith("/global/homes"):
        host_name = "cori"
    elif path_home.startswith("/home/"):
        if os.environ["SLURM_SUBMIT_HOST"] == "hpc1":
            host_name = "hpc1"
        elif os.environ["SLURM_SUBMIT_HOST"] == "hpc2":
            host_name = "hpc2"
    elif path_home.startswith("/Users/"):
        host_name = "local"
    elif path_home.startswith("/home1/"):
        host_name = "stampede"
    else:
        print("Check cluster settings")
        sys.exit()
    return host_name


def identify_vasp_environ(hpc, gamma_only):
    if hpc.lower() == "hpc1":
        os.environ[
            "VASP_PP_PATH"
        ] = "/home/ark245/programs/vasp5.4.4/pseudopotentials/pseudo54"
        if gamma_only:
            vasp_exe = "vasp_gam"
        else:
            vasp_exe = "vasp_std"
        os.environ["VASP_COMMAND"] = (
            "module load vasp/5.4.4pl2-vtst; NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo \" $NTASKS * $NNODES \" | bc`; echo \"num_cpu=\" $NCPU; srun -n $NCPU %s | tee -a op.vasp"
            % vasp_exe
        )

    elif hpc.lower() == "hpc2":
        os.environ[
            "VASP_PP_PATH"
        ] = "/home/sours/programs/vasp_PP"  # '/home/ark245/programs/pseudopotentials/pseudo54'
        if gamma_only:
            vasp_exe = "vasp_gam"
        else:
            vasp_exe = "vasp_std"
        os.environ["VASP_COMMAND"] = (
            "NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo \" $NTASKS * $NNODES \" | bc`; echo \"num_cpu=\" $NCPU; $(which mpirun) --map-by core --display-map --report-bindings --mca btl_openib_allow_ib true --mca btl_openib_if_include mlx5_0:1 --mca btl_openib_warn_nonexistent_if 0 --mca btl_openib_warn_no_device_params_found 0 --mca pml ob1 --mca btl openib,self,vader --mca mpi_cuda_support 0 -np $NCPU %s | tee -a op.vasp"
            % vasp_exe
        )
        print(os.environ["VASP_COMMAND"])

    elif hpc.lower() == "cori":
        os.environ["VASP_PP_PATH"] = "/global/homes/a/ark245/pseudopotentials/PBE54"
        if gamma_only:
            vasp_exe = "vasp_gam"
        else:
            vasp_exe = "vasp_std"
        os.environ["VASP_COMMAND"] = (
            "NTASKS=`echo $SLURM_TASKS_PER_NODE|tr '(' ' '|awk '{print $1}'`; NNODES=`scontrol show hostnames $SLURM_JOB_NODELIST|wc -l`; NCPU=`echo \" $NTASKS * $NNODES \" | bc`; echo \"num_cpu=\" $NCPU; srun -n $NCPU %s | tee -a op.vasp"
            % vasp_exe
        )

    elif hpc.lower() == "stampede":
        os.environ["VASP_PP_PATH"] = "/home1/05364/ark245/pseudopotentials/PBE54"
        if gamma_only:
            vasp_exe = "vasp_gam_vtst"
        else:
            vasp_exe = "vasp_std_vtst"
        os.environ["VASP_COMMAND"] = (
            "module load vasp/5.4.4; export OMP_NUM_THREADS=1;rm op.vasp; mpirun -np $SLURM_NTASKS %s | tee op.vasp"
            % vasp_exe
        )

    elif hpc.lower() == "local":
        os.environ["VASP_PP_PATH"] = "local_vasp_pp"
        if gamma_only:
            vasp_exe = "vasp_gam"
        else:
            vasp_exe = "vasp_std"
        os.environ["VASP_COMMAND"] = "local_%s" % vasp_exe
    else:
        print("Check cluster settings")
        sys.exit()
    vasp_pp_path = os.environ["VASP_PP_PATH"]
    vasp_command = os.environ["VASP_COMMAND"]
    return vasp_pp_path, vasp_command
