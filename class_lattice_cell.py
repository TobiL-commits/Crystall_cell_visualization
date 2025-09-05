import numpy as np


class LatticeCellCharacterization:
    def __init__(self, unit_cell_type: str = 'pc', translation: list[list[float]] = np.zeros(3)):
        unit_cell_type: str
        translation: list[list[float]]
    

class LatticeCellAtomPositions:
    
    def __init__(self, 
                lattice_vectors: list[list[float]] = np.array([[1,0,0],[0,1,0],[0,0,1]]), 
                atom_position: list[list[float]] = np.array([[0, 0, 0], [1,0,0], [1,0,1], [1,1,0], [1,1,1], 
                            [0,1,0], [0,1,1], 
                            [0,0,1]]), 
                atom_types: list[list[float]] = np.ones(8),
                atom_radii: list[list[float]] = np.ones(8) * 1/2,
                # More properties
                ) -> None:
        self.lattice_vectors = lattice_vectors
        self.atom_position = atom_position
        self.atom_types = atom_types
        self.atom_radii = atom_radii


class LatticeCellIdentifier:
    """
    Links a specific lattice cell characterization to atom positions.
    """
    def __init__(self, characterization: LatticeCellCharacterization, atom_positions: LatticeCellAtomPositions):
        self.characterization = characterization
        self.atom_positions = atom_positions


class Lattice:
    """
    Links same-celled lattice cell identifiers into a grid.
    """
    # needs to handle atoms at the same positions,
    # OR create another special class or use the old one 
    # for one celled applications
    pass