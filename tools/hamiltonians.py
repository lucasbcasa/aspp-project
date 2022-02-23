import numpy as np
import constants

def Hamiltonian(nanowire, k=None):
    n_wire  = nanowire.n_wire
        
    H_t = np.kron( np.diag( np.full((n_wire), 2 * nanowire.t) - nanowire.mu_profile ) , np.eye(2) ) 
                            # Creates an array of t's in case mu is just a float - np.diag argument must be 1d
    H_B = np.kron( np.eye( n_wire ), nanowire.B * constants.sigma_x )

    v = 1j * nanowire.tso * np.tensordot( nanowire.so_axis, constants.sigmas, axes=(0,0) ) - nanowire.t * np.eye(2)
    #v =  nanowire.tso * np.tensordot( nanowire.so_axis, constants.sigmas, axes=(0,0) ) - nanowire.t * np.eye(2)
    #v = -1j * ( (1j * nanowire.tso * constants.sigma_x) - (1j * nanowire.t * np.eye(2))) 

    
    ########
    
    try:
        tau = nanowire.tau
        vec = np.full(n_wire - 1, 1.0)
        vec[nanowire.n_sc - 1] = tau
        vec[nanowire.n_sc + nanowire.n_normal - 1] = tau
        v_grid = np.diag(vec, k=1)
    except AttributeError:
        v_grid = np.eye(n_wire, k=1 )

    H_v = np.kron(v_grid, v )
    
    
    ###########
    H_v += H_v.transpose().conjugate()
    
    if k is not None:
        #v_boundary = 1j * nanowire.tso * np.tensordot( nanowire.so_axis, constants.sigmas, axes=(0,0) ) - (nanowire.t * np.eye(2) *np.exp(1j*k))
        #H_v += np.kron(np.eye(n_wire,k=1-n_wire), v_boundary) + np.kron(np.eye(n_wire,k=n_wire-1), (v_boundary).transpose().conjugate())
        H_v += np.kron(np.eye(n_wire,k=1-n_wire), v*np.exp(1j*k)) + np.kron(np.eye(n_wire,k=n_wire-1), (v*np.exp(1j*k)).transpose().conjugate())

    return H_t + H_B + H_v

def HBdG(nanowire): 
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


# # Remaining to be implemented

# def Hamiltonian_sparse(nanowire, B): # Nanowire is taken to be in the shape of build_nanowire
#     (N,t,mu,sc,tso,sigma_axis,periodic_bc) = nanowire
#     H_t = sparse.kron(sparse.diags(2*t-mu), sparse.eye(2))
#     H_B = sparse.kron(sparse.eye(N), sparse.csr_matrix(B*pauli_matrices[0]))
#     v = (1j*tso*np.inner(sigma_axis, np.transpose(pauli_matrices, (1,2,0)))) + np.asarray([[-t, 0], [0, -t]])
#     H_v = sparse.kron(sparse.eye(N,k=1), sparse.csr_matrix(v))
#     H_v += H_v.transpose().conjugate()
#     
#     if(periodic_bc): 
#         phase_diff = 0*(nanowire[3][0]-nanowire[3][-1])
#         H_v += sparse.kron(sparse.eye(N,k=1-N), v*np.exp(1j*phase_diff/2)) + sparse.kron(sparse.eye(N,k=N-1), (v*np.exp(-1j*phase_diff/2)).transpose().conjugate())
#     
#     return H_t + H_B + H_v
# 
# def HBdG_sparse(nanowire, B):
#     (N,t,mu,sc,tso,sigma_axis,periodic_bc) = nanowire
#     H_0 = Hamiltonian_sparse(nanowire, B)
#     H_d = sparse.kron(sparse.diags(sc), sparse.csr_matrix(1j*pauli_matrices[1]))
#     
#     H = sparse.kron(sparse.csr_matrix([[1,0],[0,0]]), H_0)
#     H += sparse.kron(sparse.csr_matrix([[0,0],[0,1]]), -H_0.transpose())
#     H += sparse.kron(sparse.csr_matrix([[0,1],[0,0]]), H_d)
#     H += sparse.kron(sparse.csr_matrix([[0,0],[1,0]]), H_d.conjugate().transpose())
#     
#     return H
# 
