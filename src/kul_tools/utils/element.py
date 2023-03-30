from ase import atoms

def get_element_positions(atoms: atoms, element_1: str):
    # TODO Perform sanity check on element str
    
    positions = []
    for index, atom in enumerate(atoms):
        if atom.symbol == element_1:
            positions.append(index)

    return positions

