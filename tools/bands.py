import hamiltonians
import nanowires
import numpy as np
import solvers
import constants

# Given a set of parameters (mu, B, delta, alpha, a) - all single numbers, form a nanowire of a unit cell
# with periodic boundary conditions and wavenumber k
def unit_cell(k, **kwargs):
    n = kwargs.get('n',1)
    cell = nanowires.Nanowire(n_wire=n, **kwargs, k=k)
    return cell

# Given a set of parameters (mu, B, delta, alpha, t/a) - all single numbers
# Solve the Hamiltonian for all k values in the first BZ
def get_bands(k=constants.k_space, **kwargs):
    print(kwargs.keys())
    cell = unit_cell(-np.pi, **kwargs)
    bands = solvers.solve_k(cell, k=k)
    return bands