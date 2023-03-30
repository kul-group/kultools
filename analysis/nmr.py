#!/usr/bin/env python
import os
import warnings

from ase.io import read, write
from glob import glob
import sys
from ase.neighborlist import NeighborList, natural_cutoffs

from .vibrational_frequency import VibModes


class NMRFrequency(VibModes):
    nmr_outcar_keyword = (
        " (absolute, valence and core) "  # id's the line where relevant data is.
    )

    def __init__(self, OUTCAR, atoms=None, frequency_range=None):
        super().__init__(OUTCAR, atoms, frequency_range)

    def read(self):
        with open(self.OUTCAR) as file:
            lines = file.readlines()
        for i, line in enumerate(lines):
            if self.nmr_outcar_keyword in line:
                for atom_index in range(len(self.atoms)):
                    nmr_freq_index = (
                        i + 1 + atom_index
                    )  # i for line number in OUTCAR, +1 for skipping the line self.nmr_outcar_keyword and atom_index for different atoms in the structure
                    freq = self._get_freqs(
                        lines[nmr_freq_index]
                    )  # skipping the line self.nmr_outcar_keyword
                    if self.min_freq <= abs(freq) <= self.max_freq:
                        self.frequencies.append(freq)
                break

        if not self.frequencies:
            raise ValueError(
                "No frequency found in {0} for range [{1}, {2}]".format(
                    self.OUTCAR, self.min_freq, self.max_freq
                )
            )

    def _get_freqs(self, line, position=4) -> float:
        """
        Args:
            line (_type_): Line containing NMR frequency
            position (int, optional): Gets the 5th column containing iso_shift (incl. G=0 contribution). Defaults to 4.

        Returns:
            int: iso_shift
        """
        return float(line.split()[position])

    def write(self, path, mult=1):
        raise NotImplementedError("Write for NMR is not implemented")

    def get_freqs(self, tag, **kwargs) -> list[float]:
        tag_indicies = set(self.__get_element_position(tag))

        return [self.frequencies[index] for index in tag_indicies]

    def __get_element_position(self, tag):
        return super()._VibModes__get_element_position(tag)
