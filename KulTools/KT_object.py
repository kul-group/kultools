"""class structure(custom_params,ase_atoms)
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

.sh script package -> environments"""

from kt_structure import KT_structure
from ase.calculators.vasp import Vasp
from load_environ import identify_hpc_cluster, identify_vasp_eviron


class KT_Object(KT_structure):
    def __init__(
        self,
        gamma_only,
        structure,
        structure_type,
        calculator_name,
        calculation_type,
        calc_params={},
    ):

        super().__init__(structure, structure_type, calculator_name, calculation_type)

        self.gamma = gamma_only
        self.user_params = calc_params
        self.kt_params = {}

        if self.calc_name.lower() == "vasp":
            self.load_structure_type_params_vasp()
            self.load_calc_type_params_vasp()

        self.load_update_kt_params()

        self.init_ase_calc()

        self.hpc = identify_hpc_cluster()
        identify_vasp_eviron(self.hpc, self.gamma)

    def load_structure_type_params_vasp(self):

        if self.structure_type.lower() == "metal":
            self.default_calc_params.update({"sigma": 0.2, "ismear": 1})

        elif self.structure_type.lower() == "gas-phase":
            self.default_calc_params.update({"lreal": False})

        else:
            pass

    def load_calc_type_params_vasp(self):

        if self.calculation_type.lower() == "spe":
            self.default_calc_params.update({"nsw": 0})

        elif self.calculation_type.lower() == "opt_fine":
            self.default_calc_params.update(
                {"ibrion": 1, "potim": 0.05, "nsw": 50, "ediffg": -0.03}
            )

        elif self.calculation_type.lower() == "vib":
            self.default_calc_params.update({"ibrion": 5, "potim": 0.02, "nsw": 1})

        elif self.calculation_type.lower() == "solv":
            self.default_calc_params.update(
                {
                    "potim": 0.0,
                    "nsw": 5,
                    "lwave": True,
                    "lsol": False,
                    "prec": "Accurate",
                }
            )

        elif self.calculation_type.lower() == "solv-spe":
            self.default_calc_params.update(
                {
                    "potim": 0.0,
                    "nsw": 3,
                    "lwave": True,
                    "lsol": True,
                    "prec": "Accurate",
                }
            )

        elif self.calculation_type.lower() == "solv-spe-rho":
            self.default_calc_params.update(
                {
                    "potim": 0.0,
                    "nsw": 5,
                    "lrho": True,
                    "lwave": True,
                    "lsol": True,
                    "prec": "Accurate",
                    "lrhob": True,
                    "lrhoion": True,
                }
            )

        elif self.calculation_type.lower() == "md":
            self.default_calc_params.update(
                {"ibrion": 0, "nsw": 1000000, "tebeg": 298, "isif": 2, "smass": 0}
            )

        else:
            pass

    def load_update_kt_params(self):
        self.default_calc_params.update(self.user_params)
        self.kt_params = self.default_calc_params

    def init_ase_calc(self):
        if self.calc_name.lower() == "vasp":
            self.ase_calculator = Vasp()
        else:
            pass

        for key in self.kt_params:
            self.ase_calculator.set(key=self.kt_params[key])
