import os
import warnings

import numpy as np
from ase.atoms import Atoms
from numpy import linalg as LA
from ase.io.vasp import read_vasp_out
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


class VibModes(object):
    """
    Class for visualizing vibrational modes from VASP calculations.
    Reads OUTCAR from ibrion=5/6 calculation and writes .traj files
    to see the atomic motions corresponding to each mode.

    Args:
        OUTCAR (str): Path to OUTCAR file
        freuency_range (list[float]): Range of mode frequencies of interest, e.g.,
            [freq_min, freq_max].

    Usage:
    ```py
    OUTCAR = '<OUTCAR-dir>/OUTCAR'
    traj_file_path = '<tagged-atoms-dir>/tagged_atoms.traj'
    tag = 1
    modes = VibModes(OUTCAR, traj_file_path)
    modes.read()
    print(modes.get_freqs(tag))
    modes.write('<output-dir>') # will create the folder if doesn't exist

    ```
    """

    def __init__(self, OUTCAR, atoms: Atoms = None, frequency_range=None):
        self.set_range(frequency_range)
        self.OUTCAR = OUTCAR
        if atoms:
            self.atoms = atoms
        else:
            self.atoms = read_vasp_out(filename=self.OUTCAR)
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

        if not self.frequencies:
            warnings.warn(
                "No frequency found in {0} for range [{1}, {2}]".format(
                    self.OUTCAR, self.min_freq, self.max_freq
                )
            )

    def write(self, path, mult=1):
        # freq_range is already passed in the init function
        r = self.atoms.get_positions()
        atoms = self.atoms.copy()
        j = 0
        for freq, disp in zip(self.frequencies, self.displacements):
            if self.min_freq <= abs(freq) <= self.max_freq:
                i = "" if freq > 0 else "i"
                traj_folder = os.path.join(path, "modes")
                if not os.path.exists(traj_folder):
                    os.makedirs(traj_folder)
                traj_path = os.path.join(
                    traj_folder, f"{j:04d}_{abs(freq):.0f}{i}.traj"
                )
                traj = Trajectory(traj_path, "w")
                for x in np.linspace(0.0, 2 * np.pi, 20, endpoint=False):
                    atoms.set_positions(r + mult * np.sin(x) * disp)
                    traj.write(atoms)
                j += 1

    def get_freqs(self, tag, **kwargs) -> list[float]:
        # validate element_1 is correct
        selected_freq = []
        tag_indicies = set(self.__get_element_position(tag))

        for freq, disp in zip(self.frequencies, self.displacements):
            tagged_element_disps = []
            all_displacements = []

            for position in range(len(self.atoms)):
                delta = self.__get_displacement(disp, position)
                if position in tag_indicies:
                    tagged_element_disps.append(delta)
                else:
                    all_displacements.append(delta)

            if self.__is_prefered_displacements(
                all_displacements, tagged_element_disps
            ):
                selected_freq.append(freq)

        return selected_freq

    def __get_displacement(self, displacement_matrix, position):
        if len(displacement_matrix[position]) != 3:
            raise AttributeError
        return self.__get_l2_norm(displacement_matrix[position])

    def __get_element_position(self, tag):
        indicies = []
        for position in range(len(self.atoms)):
            if self.atoms[position].tag == tag:
                indicies.append(position)

        return indicies

    def __get_l2_norm(self, delta: list[float], direction_vector=None):
        if direction_vector:
            raise NotImplementedError
        return LA.norm(delta)

    @staticmethod
    def __is_prefered_displacements(all_displacements, tagged_element_disps):
        return sum(all_displacements) < sum(tagged_element_disps)

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
