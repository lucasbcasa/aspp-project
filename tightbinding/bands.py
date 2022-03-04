"""
Methods to calculate the band structure for a periodic nanowire.

unit_cell creates the unit cell of a system;
get_bands calls unit_cell and calculates it corresponding band structure.
"""
import nanowires
import solvers
import constants

def unit_cell(k, **kwargs):
    """
    Create a Nanowire object representing a unit cell of a periodic system.
    
    The unit cell is then returned.
    
    Parameters
    ----------
    k : float
        The k wavenumber to be applied to the unit cell's boundary.
    
    Other Parameters
    ----------------
    **kwargs : dict
        Optional parameters to be passed to the Nanowire constructor.
    """
    n = kwargs.get('n',1)
    cell = nanowires.Nanowire(n_wire=n, **kwargs, k=k)
    return cell

# Given a set of parameters (mu, B, delta, alpha, t/a) - all single numbers
# Solve the Hamiltonian for all k values in the first BZ
def get_bands(k=constants.k_space, **kwargs):
    """
    Construct and solve the Hamiltonian for a periodic system with k values in the first Brillouin Zone and parameters specified.
    
    The eigenvalues are returned.
    
    Parameters
    ----------
    k : 1-d array
        An arry with the k wavenumbers to be applied to the unit cell's boundary. Default value is saved in the constants module.
    
    Other Parameters
    ----------------
    **kwargs : dict
        Optional parameters to be passed to the Nanowire constructor in the unit_cell method.
    """
    cell = unit_cell(-np.pi, **kwargs)
    bands = solvers.solve_k(cell, k=k)
    return bands