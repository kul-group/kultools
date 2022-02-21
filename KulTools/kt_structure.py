from calc_object import default_calc
from ase import Atoms


class KT_structure(default_calc):
    def __init__(self, structure, structure_type, calculator_name, calculation_type):
        super().__init__()
        self.structure = structure
        self.structure_type = structure_type
        self.calc_name = calculator_name
        self.calculation_type = calculation_type

        self._check_calculator_name()
        self._load_default_calc_params()
        self._check_structure_type()
        self._check_calc_type()

    def _check_calculator_name(self):

        assert self.calc_name.lower() in ["vasp", "emt"], (
            "Unknown/Not implemented calculator name = %s" % self.calc_name
        )

    def _load_default_calc_params(self):

        if self.calc_name.lower() == "vasp":
            self.default_calc_params = self.load_default_vasp()

        else:
            pass

    def _check_structure_type(self):

        if self.calculation_type.lower() == "neb":
            assert isinstance(self.structure, list) and isinstance(
                self.structure[0], Atoms
            ), "Not an ASE atoms object"

        else:
            assert isinstance(self.structure, Atoms), "Not an ASE atoms object"

        assert self.structure_type.lower() in ["zeo", "mof", "metal", "gas-phase"], (
            "Unknown structure_type = %s"
            % self.structure_type  # Add insulating slab , conducting slab and gas-phase no need for metal
        )

    def _check_calc_type(self):

        assert self.calculation_type.lower() in [
            "spe",
            "opt",
            "opt_fine",
            "vib",
            "neb",
            "md",
            "solv",
            "solv-spe",
            "solv-spe-rho",
        ], (
            "Unknown calculation_type = %s" % self.calculation_type
        )
