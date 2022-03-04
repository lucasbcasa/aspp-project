
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
