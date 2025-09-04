import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

"""
# Plan 
## (a) First plot with a unit cell and scatter spheres 
## (b)
### 1. Method for sphere around a point (with x0, r, color)
### 2. Place spheres on the points. 

# TODOs (priority):

# 0. Draw any cell outline and any cell volume based on the end points

# 1. Make Wigner-Seitz cells and planes drawable

# 2. -  Place cell at x0, y0, z0, so a larger lattice can be drawn (then colors must also be changed to 
#       not random)
#    -  Be able to draw (primitive) cell with arbitrary points

# 3. -  Adjust atom shapes at the edge of the lattice when using large atoms
#    -  Allow user to decide if spheres are full or not full.
# 4. Add all other shapes
# 5. Reciprocal lattice
# 6. n-th Brillouin zones
# 7. Visualize symmetries
# 8 
"""

"""
NEED TO CHANGE
1.  the call "species" is wrong. If we want to follow the library for crystal 
    data, we have to use a different format: https://spglib.readthedocs.io/en/latest/variable.html#types
2.  Rewrite how the positions of the atoms are saved. It should be able to have initial conditions and
    so on, from which we can update initial atom's position, when we want to integrate a physics engine.
"""
# Cells
class lattice_cell :    # to recieve data, see https://spglib.readthedocs.io/en/latest/variable.html    """
    """"
    Class representing a lattice cell for crystal visualization.
    Stores basis vectors, atom positions, radii, species, and translation.
    """
    
    def __init__(cell, 
            basis = np.array([[1,0,0],[0,1,0],[0,0,1]]), 
            shape='pc', 
            atoms = np.array([[0, 0, 0], [1,0,0], [1,0,1], [1,1,0], [1,1,1], 
                            [0,1,0], [0,1,1], 
                            [0,0,1]]), ## as lattice vectors
            r=np.array([[1/2]]), 
            species = np.array([8]),
            translation = np.array([0,0,0])): ## maybe use np.shape and a multidimensional array for the atoms
        """
        Initialize the lattice cell.

        Parameters:
            basis (np.ndarray): 3x3 array of basis vectors.
            shape (str): Cell type ('pc', 'fcc', 'bcc', 'hex', etc.).
            atoms (np.ndarray): Positions of atoms in lattice coordinates.
            r (np.ndarray): Radii of atoms.
            species (np.ndarray): Number of atoms per species.
            translation (np.ndarray): Translation vector for the cell.
        """        
        cell.shape = shape
        cell.basis = basis
        cell.abs_a , cell.abs_b, cell.abs_c = np.array([[np.sqrt(np.sum(basis[0]**2))], [np.sqrt(np.sum(basis[1]**2))], [np.sqrt(np.sum(basis[2]**2))]])
        cell.atoms = atoms
        cell.r = r
        cell.translation = translation
        if len(r) != len(species) :
            print("Radii do not correspond with different species. It is assumed they are the size.")
            cell.r = r
            for e in species[1:] : 
                cell.r = np.append(cell.r, r[0])

        
        if np.sum(species) == len(atoms[:]):
            cell.species = species
        else : 
            print("Invalid species count. It will be assumed, all atoms are of the same species.")
            cell.species = np.array([len(atoms[:])])

    def draw_basis(cell, ax, order=0):
        """
        Draw the basis vectors of the cell.

        Parameters:
            ax: Matplotlib 3D axis.
            order (int): Z-order for drawing.
        """
        a,b,c = cell.basis
        origin = cell.translation
        ax.quiver(*origin, a[0], 
                a[1], 
                a[2], 
                color='red', label='a',arrow_length_ratio=0.1, zorder = 5)
        ax.quiver(*origin, b[0],  b[1], b[2],  color='green', label='b',arrow_length_ratio=0.1,zorder = 5)
        ax.quiver(*origin, c[0],c[1], c[2],  color='blue', label='c',arrow_length_ratio=0.1,zorder = 0)
        
        
    # Draw all edges and return to origin
    def cell_edges(cell):
        """
        Compute the edges of the cell for plotting.

        Returns:
            Arrays of x, y, z coordinates for the edges.
            For 'fcc' and 'bcc', also returns additional edge arrays.
        """
        a, b, c = cell.basis

        # Define the vertices of the shape based on the basis vectors
        if cell.shape == 'hex': 
            ## insert hex edges here
            pass 
        else : 
            vertices = np.array([
                [0, 0, 0],
                a,
                a + b,
                b,
                c,
                c + a,
                c + a + b,
                c + b
            ])

            # Define the order to connect the vertices to trace the edges
            edge_indices = [(0, 1), (2, 3), (0, 4), (5, 6), (7, 4), (0, 1), (5, 6), (2, 3), (7,4), (0,0)]

            # Initialize lists to store the x, y, and z coordinates of the traced edges
            x = []
            y = []
            z = []

            # Trace the edges by connecting the vertices in the specified order
            for i, j in edge_indices:
                x.extend([vertices[i, 0], vertices[j, 0]])
                y.extend([vertices[i, 1], vertices[j, 1]])
                z.extend([vertices[i, 2], vertices[j, 2]])

            #conditionals if cell type is specified further   
            if cell.shape == 'pc':
                return np.array(x),np.array(y),np.array(z) 
            else:
                    # Initialize lists to store the x, y, and z coordinates of the traced edges
                l = []
                m = []
                n = []
                if cell.shape == 'fcc': 
                    edge_indices = [(0,2), (7,5),(0,3),(1,6),(4,3),(6,2),(5,1),(4,0),(5,0)]       
                elif cell.shape == 'bcc':
                    edge_indices = [(0,6), (5,3), (7,1), (2,4), (0,0)]
                for i, j in edge_indices:
                    l.extend([vertices[i, 0], vertices[j, 0]])
                    m.extend([vertices[i, 1], vertices[j, 1]])
                    n.extend([vertices[i, 2], vertices[j, 2]])
                return np.array(x),np.array(y),np.array(z), np.array(l),np.array(m),np.array(n)

    def draw_cell(cell, ax, order = 1000,            # order set high to be drawn over the spheres 
                opacity = 1,color = 'k'):  
        """
        Draw the outline of the cell.

        Parameters:
            ax: Matplotlib 3D axis.
            order (int): Z-order for drawing.
            opacity (float): Opacity of the lines.
            color (str): Color of the cell edges.
        """
        out = None
        
        if (cell.shape == 'fcc') | (cell.shape == 'bcc'):

            x,y,z, l,m,n = cell.cell_edges()

            out = ax.plot(x + cell.translation[0] * np.ones(x.shape) ,y + cell.translation[1] * np.ones(y.shape),z + cell.translation[2] * np.ones(z.shape),
                        'k',zorder=order,alpha=opacity), ax.plot(l + cell.translation[0] * np.ones(l.shape), m + cell.translation[1] * np.ones(m.shape), n + cell.translation[2] * np.ones(n.shape),
                        'k:', zorder=order, alpha=opacity,color = color)
        else:
            x,y,z = cell.cell_edges()
            out = ax.plot(x + cell.translation[0] * np.ones(x.shape) ,y + cell.translation[1] * np.ones(y.shape),z + cell.translation[2] * np.ones(z.shape),
                        'k-',zorder=order, alpha=opacity,color = color)
        return out
    
    def draw_sphere(cell, ax, order,
            r, pos = np.array([[1/2,1/2,1/2]]), dir = [[1,0,0],[0,1,0],[0,0,1]],     # pos in units of 1/a ## defining the direction of the sphere (## bulge of partial ball in the direction of the vector dir)
            frac = 'full', size='big', opacity = 1, color = 'g'):                                  # frac is the fraction of the ball ('full', 'half' , 'other' == the shape in the vectors)
        # We make one full ball made up of points 
        # we then alter the shape of the intervals based on its frac and dir and generate boundary durvaces on the other side of the ball
        """
        Draw a sphere (atom) at a given position.

        Parameters:
            ax: Matplotlib 3D axis.
            order (int): Z-order for drawing.
            r (float): Radius of the sphere.
            pos (np.ndarray): Position of the sphere (fractional coordinates).
            dir (list): Directions for partial spheres (not implemented).
            frac (str): Fraction of the sphere to draw ('full', 'half', 'other').
            size (str): 'big' or 'small' for scaling.
            opacity (float): Opacity of the sphere.
            color: Color of the sphere.
        """
        out = None
        n_of_tiles = 100

        if size =='small':
            r *= 0.2
        elif size == 'big':
            pass

        # generate values for a sphere
        u = np.linspace(0, 2 * np.pi, n_of_tiles)    # phi, for an eigth we need 2pi/4
        v = np.linspace(0, np.pi, n_of_tiles)        # theta, for an eigth, we need pi/2=2pi/4

        x = r * np.outer(np.cos(u), np.sin(v)) + pos[0] * cell.abs_a 
        y = r * np.outer(np.sin(u), np.sin(v)) + pos[1] * cell.abs_b 
        z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + pos[2] * cell.abs_c 

        x += cell.translation[0] * np.ones(x.shape)
        y += cell.translation[1] * np.ones(y.shape)
        z += cell.translation[2] * np.ones(z.shape)

        # ToDo: Cut and orient the full sphere according to dir 
        if frac == 'half':      #Idea: see, if some values of the combined vectors are in some volume, then delete them. Afterwards, append all the area values for the closing surfaced to the xyz-arrays
            pass
            
        elif frac == 'other': #schauen ob punkte innerhalb der vektoren sind. 
            pass
        elif frac == 'full':
            pass

        #print the surface of the xyz-values.
        out = ax.plot_surface(x, y, z, zorder=order, color=color,alpha = opacity)   
        out = out, ax.set_aspect('equal')
        return out
    
    def draw_atoms(cell, ax,
            order=0, 
            size = 'big', 
            opacity = 1):  ### The Other shapes are not included yet
        """
        Draw all atoms in the cell.

        Parameters:
            ax: Matplotlib 3D axis.
            order (int): Z-order for drawing.
            size (str): 'big' or 'small' for scaling.
            opacity (float): Opacity of the atoms.
        """
        niterator = 0 
        riterator = 0
        if size == 'small':
            cell.r *= 0.5
        
        for n in cell.species:

            color = np.random.rand(3)
            
            jiterator = 0
            for m in cell.atoms[niterator:n+niterator]:           #iterate from the nth element to the last element of it's species
                cell.draw_sphere(ax, order+(n+1)*5, r = cell.r[riterator], pos = cell.atoms[niterator + jiterator], color = color.tolist(), size=size, opacity=opacity)
                jiterator +=1
            niterator += n
            riterator +=1
        if size == 'small' : 
            cell.r *= 1/0.5