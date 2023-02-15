import ase
from ase.geometry.analysis import Analysis

from .element import get_element_positions

def set_tags(atoms: ase.Atoms, tag: int, element_1: str, element_2: str, do_sanity_check: bool=True):
    """Set tags in `ase.atoms` according to bonds

    Args:
        atoms (ase.Atoms): ase atoms object of the molecule/structure
        tag (int): tag to set to the specified bonded atoms
        element_1 (str): element symbol eg N, C, O, Ti
        element_2 (str): element symbol eg N, C, O, Ti
        do_sanity_check (bool, optional): whether to perform sanity check?. Defaults to True.
    """
    def tag_sanity_check(atoms):
        if len(set(atoms.get_tags())) > 1:
            raise ValueError("Expected no initial tags set")

    if do_sanity_check:
        tag_sanity_check(atoms)

    element_1_positions = get_element_positions(atoms, element_1)

    analysis = Analysis(atoms)
    all_bonds = analysis.all_bonds[0]

    for element_1_position in element_1_positions:
        for bonded_atom_position in all_bonds[element_1_position]:
            if atoms[bonded_atom_position].symbol == element_2:
                atoms[element_1_position].tag = tag
                atoms[bonded_atom_position].tag = tag

