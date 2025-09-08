import numpy as np


""" (PROBABLY) TO BE PUT INTO OTHER CLASSES """

def unit_vector(vector: list[float]) -> list[float]:
    return vector / np.linalg.norm(vector)

def law_of_cosines(vector_a: list[float], vector_b: list[float]) -> float:
    vector_a, vector_b = unit_vector(vector_a), unit_vector(vector_b)
    return np.arccos(np.clip(np.dot(vector_a, vector_b), -1.0, 1.0))



class LatticeCellAtomPositions:
    def __init__(
        self,
        lattice_vectors: list[list[float]] = np.array([
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ]),
        atom_position: list[list[float]] = np.array([
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [1, 0, 1],
            [1, 1, 1],
            [0, 1, 0],
            [0, 1, 1],
            [0, 0, 1]
        ]),
        atom_types: list[list[float]] = np.ones(8),
        atom_radii: list[list[float]] = np.ones(8) * 1/2,
        # ?? More properties
    ):
        self.lattice_vectors = lattice_vectors
        self.atom_position = atom_position
        self.atom_types = atom_types
        self.atom_radii = atom_radii
        
    def __repr__(self):
        return (
            f'LatticeCellAtomPositions('
            f'{self.lattice_vectors}, '
            f'{self.atom_position}, '
            f'{self.atom_types}, '
            f'{self.atom_radii})'
        ) 
    

class CalculateLatticeCrystalSystem:
    """
    Helper class to figure out what crystal structure is given.
    Important to later chose a dedicated tracing of lines to visualize the
    cell borders.
    """
    # !! Check if names are correct 
    def is_cubic_cell(atoms) -> bool:
        pass
    
    def is_tetragonal_cell(atoms) -> bool:
        pass
    
    def is_orthorhombic_cell(atoms) -> bool:
        pass
    
    def is_hexagonal_cell(atoms) -> bool:
        pass
    
    def is_rhombohedral_cell(atoms) -> bool: # also counts as hex
        pass
    
    def is_monoclinic_cell(atoms) -> bool:
        pass
    
    def is_triclinic_cell(atoms) -> bool:
        pass
    
    def calc_unit_cell_type(self) -> str: 
        """
        See Gross, Marx: Festkörperphysik De Gruytre 2022; pages 16 - 19.
        """
        if self.is_cubic_cell(self.atoms):
            return 'pc'
        ...
    

class LatticeCellCharacterization:
    def __init__(self, 
                atoms: LatticeCellAtomPositions, 
                unit_cell_type: str = 'pc', 
                translation: list[list[float]] = np.zeros(3)
                ):
        self.atoms = atoms      # !! I will probably only need the lattice_vectors not the whole class
        self.unit_cell_type = unit_cell_type
        self.translation = translation


class LatticeCellIdentifier:
    """
    Links a specific lattice cell characterization to atom positions and gives 
    cells an identity. 
    Useful also when highlighting specific parts of the Plots or for example 
    drawing a plane for a cell that is different from the origin. 
    
    !! Need to later make a funtion using plotly, to display numbers 
    and also translation vectors of the different cells within the plot.
    """
    def __init__(self, 
                characterization: LatticeCellCharacterization, 
                atom_positions: LatticeCellAtomPositions
                ):
        self.characterization = characterization
        self.atom_positions = atom_positions


class Lattice:
    """
    Links same-celled lattice cell identifiers into a grid and passing them with an Identity into 
    """
    # ?? needs to handle atoms at the same positions,
    # OR create another special class or use the old one 
    # for one celled applications
    def __repr__(self):
        pass                # ?? write function, that passes an ordered description of the cells



def main():
    pass

if __name__ == "__main__":
    main()