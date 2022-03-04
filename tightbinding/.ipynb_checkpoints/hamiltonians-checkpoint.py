"""
Construct the Hamiltonians for a given system.

Hamiltonian builds the normal state Hamiltonian.
HBdG builds the BdG Hamiltonian.
"""
import numpy as np
import constants

def Hamiltonian(nanowire, k=None):
    """
    Builds the normal state Hamiltonian for a given Nanowire object.
    
    The Nanowire object can be a Nanowire subclass such as a SNS_Junction object.
    
    Parameters
    ----------
    nanowire : nanowires.Nanowire object
        The Nanowire object representing the system. All physical parameters are taken from the Nanowire object's attributes.
    
    Other Parameters
    ----------------
    k : float, optional
        The wavenumber of a periodic system. If `k` is given it assumes the system is periodic and treats `nanowire` as the unit cell of the system.
    
    Returns
    -------
    H : np.ndarray
        The dense-matrix representation of the normal state Hamiltonian. The ndarray's shape is (2`N`,2`N`), where `N` is the length of the nanowire provided, `nanowire.n_wire`. The factor of two comes from the spin DOF.
    """
    n_wire  = nanowire.n_wire
        
    H_t = np.kron( np.diag( np.full((n_wire), 2 * nanowire.t) - nanowire.mu_profile ) , np.eye(2) ) 
                            # Creates an array of t's in case mu is just a float - np.diag argument must be 1d
    H_B = np.kron( np.eye( n_wire ), nanowire.B * constants.sigma_x )

    v = 1j * nanowire.tso * np.tensordot( nanowire.so_axis, constants.sigmas, axes=(0,0) ) - nanowire.t * np.eye(2)
    
    try:
        tau = nanowire.tau
        vec = np.full(n_wire - 1, 1.0)
        vec[nanowire.n_sc - 1] = tau
        vec[nanowire.n_sc + nanowire.n_normal - 1] = tau
        v_grid = np.diag(vec, k=1)
    except AttributeError:
        v_grid = np.eye(n_wire, k=1 )

    H_v = np.kron(v_grid, v )
    
    H_v += H_v.transpose().conjugate()
    
    if k is not None:
        H_v += np.kron(np.eye(n_wire,k=1-n_wire), v*np.exp(1j*k)) + np.kron(np.eye(n_wire,k=n_wire-1), (v*np.exp(1j*k)).transpose().conjugate())

    return H_t + H_B + H_v

def HBdG(nanowire):
    """
    Builds the Bogoliubov-de Gennes Hamiltonian for a given Nanowire object.
    
    The Nanowire object can be a Nanowire subclass such as a SNS_Junction object.
    
    Parameters
    ----------
    nanowire : nanowires.Nanowire object
        The Nanowire object representing the system. All physical parameters are taken from the Nanowire object's attributes.
    
    Returns
    -------
    HBdG : np.ndarray
        The dense-matrix representation of the Bogoliubov-de Gennes Hamiltonian. The ndarray's shape is (4`N`,4`N`), where `N` is the length of the nanowire provided, `nanowire.n_wire`. A factor of two comes from the spin DOF and the other from the particle-hole DOF.
    """
    delta_profile = nanowire.delta_profile

    H_0 = Hamiltonian(nanowire, k=nanowire.k)
    H_d = np.kron(np.diag(delta_profile), 1j * constants.sigma_y)
    #H_d = np.kron(np.diag(delta_profile), -1j * (1j * constants.sigma_x))
    
    H =  np.kron( np.asarray([[1,0],[0,0]]), H_0 )
    
    if nanowire.k is None:
        H_hole = -H_0.transpose() # Instead of this I need to calculate H_0 for -k and use that
    else:
        H_hole = -Hamiltonian(nanowire, k=-nanowire.k).transpose()
    H += np.kron( np.asarray([[0,0],[0,1]]), H_hole )
    H += np.kron( np.asarray([[0,1],[0,0]]), H_d )
    H += np.kron( np.asarray([[0,0],[1,0]]), H_d.conjugate().transpose() )
    
    return H