from copy import copy


class run_KulTools(KT_Object):

    def __init__(self,gamma_only,structure,structure_type,calculator_name,calculation_type,calc_params={}):

        super().__init__(gamma_only,structure,structure_type,calculator_name,calculation_type,calc_params)

        self.run_calculation()

            
    def run_calculation(self):

        if self.calculation_type.lower() == 'opt':
            self.run_opt()
        elif
        elif
        else

    def run_opt(self):
        ase_atoms = copy.deepcopy(self.structure)

        dir_name = 'opt'

        #Set working directory generations
        self.ase_calculator.set('directory':os.getwd())

        ase_atoms.set_calculator(self.ase_calculator)
        energy = atoms.get_potential_energy()

        
        new_atoms = self.run_dft(atoms,dir_name)
        return new_atoms

