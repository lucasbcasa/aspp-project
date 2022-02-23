import numpy as np

# Builds a probe connected to the middle of the nanowire:
def make_probe_hamiltonian(N_wire, N_probe, t, mu, B, tso, periodic_bc=False, probe_pos=0.5, strength=1/20, sigma_axis=[0,0,1]): 
    if((len(t[0])!=N_wire)|(len(mu[0])!=N_wire)|(len(B[0])!=N_wire)|(len(tso[0])!=N_wire)):
        print("Nanowire parameters not compatible")
        return
    
    if((len(t[1])!=N_probe)|(len(mu[1])!=N_probe)|(len(B[1])!=N_probe)|(len(tso[1])!=N_probe)):
        print("Probe parameters not compatible")
        return
    
    H_wire = make_hamiltonian(N_wire, t[0], mu[0], B[0], tso[0], periodic_bc, )
    
    H_probe = make_hamiltonian(N_probe, t[1], mu[1], B[1], tso[1], False, sigma_axis)
    
    # Connecting probe and wire:
    #Does tso connect probe and wire?
    itdoes = False
    v = np.asarray([[-t[1]*strength, tso[1]*strength*itdoes], [-tso[1]*strength*itdoes, -t[1]*strength]])
    H_v = np.zeros((2*N_wire,2*N_probe))
    probe_index = int(N_wire*probe_pos)
    H_v[2*probe_index][0] = v[0][0][0]
    H_v[2*probe_index][0+1] = v[0][1][0]
    H_v[2*probe_index+1][0] = v[1][0][0]
    H_v[2*probe_index+1][0+1] = v[1][1][0]
    
    H1 = np.concatenate((H_wire, H_v), axis=1)
    H2 = np.concatenate((H_v.conj().T, H_probe), axis=1)
    H = np.concatenate((H1,H2))
    
    return H

# Builds the complete BdG Hamiltonian:
def make_ham_hole_sym_probe(N_wire, N_probe, t, mu, B, tso, sc, periodic_bc=False, probe_pos=0.5, strength=1/20, sigma_axis=[0,0,1]): 
    
    
    
    if((len(t[0])!=N_wire)|(len(mu[0])!=N_wire)|(len(B[0])!=N_wire)|(len(tso[0])!=N_wire)):
        print("Nanowire parameters not compatible")
        return
    
    if((len(t[1])!=N_probe)|(len(mu[1])!=N_probe)|(len(B[1])!=N_probe)|(len(tso[1])!=N_probe)):
        print("Probe parameters not compatible")
        return
    
    H0 = make_probe_hamiltonian(N_wire, N_probe, t, mu, B, tso, periodic_bc, probe_pos, strength, sigma_axis)
    
    N = N_wire + N_probe
    profile = np.concatenate((sc[0],sc[1]))
    d = np.asarray([[np.zeros(N), profile],[-profile, np.zeros(N)]]) 
    # @Could maybe use kronecker products with this matrix to build H
    H = np.zeros((4*N,4*N), dtype=complex)
    
    for i in range(2*N): #adding h
        for j in range(2*N):
            
            H[i][j] = H0[i][j]
            H[i+2*N][j+2*N] = -np.conjugate(H0[i][j])

    for i in range(N) :
        H[2*i][2*N+2*i+1] = d[0][1][i]
        H[2*i+1][2*N+2*i] = d[1][0][i] 
    for i in range(N) :
        H[2*N+2*i][2*i+1] = np.conjugate(d[1][0][i])
        H[2*N+2*i+1][2*i] = np.conjugate(d[0][1][i])
    return H

def solve_field_probe(N_wire, N_probe, t, mu, Bc, tso, sc, b=np.linspace(0,2,200), periodic_bc=False, probe_pos=0.5,strength=1/20):
    steps = len(b)
    spec = []
    start = time.time()
    for i in range(steps):
        energies = np.linalg.eigvalsh(make_ham_hole_sym_probe(N_wire, N_probe, t, mu, b[i]*Bc, tso, sc, periodic_bc, probe_pos,strength))
        spec.append(energies)
        if (10*i in np.linspace(0,9*steps,10, dtype=int)):
            print(round(100*i/steps,1),"%", end ='|')
    end = time.time()
    print("Time elapsed when solving for", N_wire+N_probe, "sites: ", (end - start)/60, "min")
    return np.asarray(spec)


def plot_wavefunction_probe(N_wire, N_probe, t, mu, Bc, tso, sc, b, states=[0], periodic_bc=False, probe_pos=0.5,strength=1/20, sigma_axis=[0,0,1]): # b is an array of values
    wav =[] 
    spec = []
    for i in range(len(b)):
        eigs = np.linalg.eigh(make_ham_hole_sym_probe(N_wire, N_probe, t, mu, b[i]*Bc, tso, sc, periodic_bc, probe_pos,strength, sigma_axis))
        spec.append(eigs[0])
        wav.append(np.transpose(eigs[1]))
    color = np.transpose([np.linspace(0,1,len(states)),np.zeros(len(states)),np.linspace(1,0,len(states))])
    plt.figure(num=None, figsize=(6,4), dpi=200, facecolor='w', edgecolor='k')
    for i in range(len(b)):
        ax = plt.subplot(1,len(b),i+1,projection='3d')
        ax.set_yticks([0, int(N_probe/2), int(N_probe)])
        ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
        plt.title('$B=${}$B_c$'.format(b[i]))
        for j in range(len(states)):
            amplitude = wavefunction(wav[i][(2*(N_wire+N_probe)+states[j])], N_probe)
            amplitude[1] = np.insert(amplitude[1], 0, amplitude[0][int(probe_pos*N_wire)])
            print(spec[i][2*(N_wire+N_probe)+states[j]])
            #print(amplitude[0], amplitude[1])
            #print(len(amplitude[0]), len(amplitude[1]))
            plt.plot(np.linspace(0,N_wire-1,N_wire), amplitude[0], zs=0, zdir='y', c=color[j], lw=1)
            plt.plot(np.linspace(0,N_probe,N_probe+1), amplitude[1], zs=int(probe_pos*N_wire), zdir='x', c=color[j], lw=1)
        plt.ylabel("Probe wire [a]")
        plt.xlabel("Main wire [a]")
        ax.set_zlabel(r'$\psi^2$')
        ax.view_init(elev=20, azim=135)
        ax.dist = 13
        ax.invert_yaxis()
        ax.set_xticks([int(N_wire), int(N_wire/2), 0])
        plt.legend((blue,red),("Lowest State","1st Excited State"))
    return