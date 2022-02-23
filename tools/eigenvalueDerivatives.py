import numpy as np
import constants

# Works for any nanowire geometry
def H_sc(nanowire):
    delta_profile = nanowire.delta_profile
    
    H_d = np.kron(np.diag(delta_profile), 1j * constants.sigma_y)
    
    H = np.kron( np.asarray([[0,1],[0,0]]), H_d )
    H += np.kron( np.asarray([[0,0],[1,0]]), H_d.conjugate().transpose() )
    return H

# Works for SNS junctions only
def H_prime(junction):
    # Delta(phi) = Delta * [e^(-i phi), 0, 1] -> dDelta/dphi = Delta * [-ie^(-i phi), 0, 0]
    delta_profile_prime = np.asarray([(-1j * d) if i<junction.n_left else 0 for i, d in enumerate(junction.delta_profile) ])
    
    H_d = np.kron(np.diag(delta_profile_prime), 1j * constants.sigma_y)
    
    H = np.kron( np.asarray([[0,1],[0,0]]), H_d )
    H += np.kron( np.asarray([[0,0],[1,0]]), H_d.conjugate().transpose() )
    return H

def H_double_prime(junction):
    
    # Delta(phi) = Delta * [e^(-i phi), 0, 1] -> dDelta/dphi = Delta * [-ie^(-i phi), 0, 0]
    delta_profile_double_prime = np.asarray([(-d) if i<junction.n_left else 0 for i, d in enumerate(junction.delta_profile) ])
    
    H_d = np.kron(np.diag(delta_profile_double_prime), 1j * constants.sigma_y)
    
    H = np.kron( np.asarray([[0,1],[0,0]]), H_d )
    H += np.kron( np.asarray([[0,0],[1,0]]), H_d.conjugate().transpose() )
    return H
   

# def dE_dphi(H, H_prime):
    
#     _ , eigenstates = np.linalg.eigh(H)
    
#     dE_dphi = np.diag( (eigenstates.conjugate().T) @ H_prime @ (eigenstates) )

#     return dE_dphi

# Gets array of sorted eigenvalues, returns iterator of tuples for degenerate subspaces indices and counts [(i0, c0), (i1,c11), ...]
def degenerate_subspaces(eigenvalues):
    
    _, idx_start, count = np.unique(eigenvalues, return_counts=True, return_index=True)

    idx_start = idx_start[count > 1]
    count = count[count > 1]
    
    return zip(idx_start, count)
    

def fix_degeneracy(eigenvalues, eigenvectors, E_prime, U, debug=False):
    # Handling degeneracies:
    
    
    
    if debug: print('Solving degeneracy issues...')
    
    # We want to diagonalize H_prime in the degenerate subspaces of H
    for index, count in degenerate_subspaces(eigenvalues):
        
        if debug: print(index, count)
        # Build and diagonalize E_prime in the given degenerate subspace of H
        H_prime_sub = np.copy(E_prime[index:index+count, index:index+count])
        _, eigenvecs_sub = np.linalg.eigh(H_prime_sub)
        
        # Calculate the unitary matrices in full space
        U_sub = np.eye(len(U))
        U_sub[index:index+count, index:index+count] = eigenvecs_sub
        if debug: print(U_sub == np.eye(len(U)))
        U = U @ U_sub
        
    return U
    
def first_order(eigenvalues, eigenvectors, H_prime, degenerate=False, debug=False):
    
    E_prime = eigenvectors.T.conjugate() @ H_prime @ eigenvectors
        
    U = np.eye(len(eigenvectors))
    
    if degenerate: U = fix_degeneracy(eigenvalues, eigenvectors, E_prime, U, debug)
    
    eigenvalues_prime = np.diag( U.T.conjugate() @ E_prime @ U )
    
    if debug:
        M = U.T.conjugate() @ E_prime @ U
        print(np.count_nonzero(M - np.diag(np.diagonal(M))) is 0)
        print(U)
    
    return eigenvalues_prime

def second_order(eigenvalues, eigenvectors, H_prime, H_double_prime, **kwargs):
    
    E_prime = eigenvectors.T.conjugate() @ H_prime @ eigenvectors
    eigenvalues_prime = np.diag(E_prime)
    
    E_mat = np.repeat(eigenvalues[None,...],len(eigenvalues),axis=0)
    E_mat -= E_mat.T
    
    A = np.divide(E_prime, E_mat, out=np.zeros_like(E_prime), where=E_mat!=0)
    
    # B is the difference between the matrix product and element-wise product between E_prime and A
    B = 2 * np.diag( E_prime @ A )
    
    #E_double_prime = (eigenvectors.T.conjugate() @ H_double_prime @ eigenvectors) + np.diag(B)

    eigenvalues_double_prime = np.diag(eigenvectors.T.conjugate() @ H_double_prime @ eigenvectors) + B
    
    return eigenvalues_double_prime