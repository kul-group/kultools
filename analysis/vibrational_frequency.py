import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from ase.io import iread, write
from ase.geometry.analysis import Analysis
from ase.io.trajectory import Trajectory

"""
TODO:
- Update documentation
- optimise read
    - use numpy array
    - maybe use regex for freq search
    - stop enumerate when no all freq read
- Check logic for get_freq(element)
- Add sanity checks
"""

class VibModes:
    """
    Class for visualizing vibrational modes from VASP calculations.
    Reads OUTCAR from ibrion=5/6 calculation and writes .traj files
    to see the atomic motions corresponding to each mode.

    Args:
        OUTCAR (str): Path to OUTCAR file
        freuency_range (list[float]): Range of mode frequencies of interest, e.g.,
            [freq_min, freq_max].
    """

    def __init__(self, OUTCAR, frequency_range=None):
        self.set_range(frequency_range)
        self.OUTCAR = OUTCAR
        self.atoms = next(iread(OUTCAR, index=0))
        self.n_atoms = len(self.atoms)
        self.chemical_symbols = self.atoms.get_chemical_symbols()
        self.frequencies = []
        self.displacements = []

    def set_range(self, frequency_range):
        if frequency_range:
            self.min_freq, self.max_freq = frequency_range
        else:
            self.min_freq = -np.inf
            self.max_freq = np.inf

    def read(self, read_displacements=True):
        with open(self.OUTCAR) as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
            if "THz" in line:
                freq = self._get_freq(line)
                if self.min_freq <= abs(freq) <= self.max_freq:
                    self.frequencies.append(freq)
                    if read_displacements:
                        disp = self._get_disp(lines[i + 2 : i + self.n_atoms + 2])
                        self.displacements.append(disp)

    def write(self, mult=1):
        # freq_range is already passed in the init function
        r = self.atoms.get_positions()
        atoms = self.atoms.copy()
        j = 0
        for freq, disp in zip(self.frequencies, self.displacements):
            if self.min_freq <= abs(freq) <= self.max_freq:
                i = "" if freq > 0 else "i"
                traj = Trajectory(f"modes/{j:04d}_{abs(freq):.0f}{i}.traj", "w")
                for x in np.linspace(0.0, 2 * np.pi, 20, endpoint=False):
                    atoms.set_positions(r + mult * np.sin(x) * disp)
                    traj.write(atoms)
                j += 1

    def get_freqs(self, element_1: str, element_2=None, **kwargs) -> list[int]:

        # validate element_1 is correct

        prim_elmt_posns = set(self.__get_element_position(element_1))

        element_1_freqs = []

        for freq, disp in zip(self.frequencies, self.displacements):
            primary_elmt_disps = []
            all_displacements = []

            for position in range(len(disp)):
                delta = self.__get_displacement(disp, position)
                l2_norm = self.__get_l2_norm(delta)
                all_displacements.append(l2_norm)
                if position in prim_elmt_posns:
                    primary_elmt_disps.append(l2_norm)

            if self.__is_prefered_displacements(all_displacements, primary_elmt_disps):
                element_1_freqs.append(freq)

        return element_1_freqs

    def __get_displacement(self, displacement_matrix, position):
        if len(displacement_matrix[position]) != 3:
            raise AttributeError
        return displacement_matrix[position]

    def __get_element_position(self, element):
        element_positions = []
        chemical_symbols = self.chemical_symbols
        for i in range(len(chemical_symbols)):
            if chemical_symbols[i] == element:
                element_positions.append(i)

        return element_positions

    def __get_bonds(self, primary_element_position, secondary_element_positions):

        preferred_bonds = []
        atom = self.atoms[primary_element_position]
        analysis = Analysis(atom.atoms)
        all_bonds = analysis.all_bonds[0]
        primary_atom_bonds = all_bonds[primary_element_position]
        for primary_atom_bond in primary_atom_bonds:
            if primary_atom_bond in secondary_element_positions:
                preferred_bonds.append(primary_atom_bond)

        return preferred_bonds

    def __get_l2_norm(self, delta: list[float], direction_vector=None):
        if direction_vector:
            raise NotImplementedError
        return LA.norm(delta)

    @staticmethod
    def __is_prefered_displacements(all_displacements, primary_elmt_disps):
        return max(all_displacements) in primary_elmt_disps

    @staticmethod
    def _get_freq(line):
        i_freq = 6 if "f/i" in line else 7
        imag = -1 if "f/i" in line else 1
        return float(line.split()[i_freq]) * imag

    @staticmethod
    def _get_disp(lines):
        def _str_to_float(line):
            return list(map(float, line.split()[-3:]))

        return np.array([_str_to_float(l) for l in lines])
