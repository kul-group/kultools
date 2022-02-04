class structure(custom_params,ase_atoms)
    assert param has been Changed
    ISMEAR SIGMA etc
    structure.get_calc_object (VASP,QE)
    structure.set_defaults

class calculation_type(custom_params):
    self.params = {}
    self.params["IBRION"] = self.custom
    contains default params for each calc
    if calculator_name - 'VASP':
        kt_params = ""
    else assert
assert param has been Changed
IBRION,NSW,EDIFFG,POTIM


from kt.calculation_type import calc_params...

.sh script package -> environments

class KT_Object:

    def __init__(self,structure,calc_params): 

        kt_params = {IBRION,NSW,EDIFF,etc}
        self.params = kt_params
        self.structure = structure
        self.calculation_type = calculation_type

    def load_params(self):

        self.hpc = self.identify_hpc_cluster()
        print('KT: HPC= %s' %self.hpc)
        self.gamma_only = gamma_only
        print('KT: VASP_GAMMA= %s' %self.gamma_only)
        
        self.structure_type = structure_type
        self.calculation_type = calculation_type
        self.main_dir = os.getcwd()
        self.structure = structure 

        self.identify_vasp_eviron()
        print('KT: VASP_PP_PATH= %s' %self.vasp_pp_path)
        print('KT: VASP_COMMAND= %s' %self.vasp_command)

        self.assign_default_calculator()