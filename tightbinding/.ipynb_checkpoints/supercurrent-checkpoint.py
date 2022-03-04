# Not adapted for temperature

import numpy as np

def derivative_of_curve(f, x=[], periodic = True):
    if(not len(x)): x=np.linspace(0,2*np.pi,len(f),endpoint=True)
    f_prime = np.zeros(len(f))
    for i in range(-1,len(f)-1):
        f_prime[i] = (f[i+1] - f[i-1])/((x[i+1]-x[i-1]))
    if periodic:
        l = x[-1]-x[0]+x[1]
        f_prime[0] = (f[1] - f[-1])/(x[1]-(x[-1]-l))
        f_prime[-1] = (f[0] - f[-2])/(x[0]-(x[-2]-l))
        
    return f_prime

def intercalate(a,b):
    if(not(len(a)==len(b))): 
        if(not(len(a)==(len(b)+1))):
            print("Failure")
            return
    #c = np.empty((a.size + b.size,), dtype=a.dtype)
    c = [[] for i in range(len(a)+len(b))]
    c[0::2] = a
    c[1::2] = b
    return c

def half_point(x):
    return (x[len(x)//2]+x[len(x)//2-1])/2

# Takes BdG spectrum as function of phase diff. between outer S reg.
# When calculating in the context of SNSNS, this will give double the correct value of current
def supercurrent(spectrum, cutoff=None, split=False): 
    
    # Filtering out negative part of the spectrum:
    if split: 
        spectrum = np.split(np.copy(spectrum),2,axis=1)[1]
    
    (N,M) = spectrum.shape
    if(cutoff): M = cutoff
    part = np.zeros((N,M))
    for i in range(N-1):
        for j in range(M):
            part[i][j] += -(spectrum[i+1][j]-spectrum[i-1][j])*N/(4*np.pi)
    for j in range(M):
        part[N-1][j] += -(spectrum[0][j]-spectrum[N-2][j])*N/(4*np.pi)
    total = np.sum(part, axis=1)
    return [np.transpose(part,axes=[1,0]), total]

def supercurrent_prime_spec(spectrum, cutoff=None, split=False): # three-point second derivative of spectrum
    
    # Filtering out negative part of the spectrum:
    if split: 
        spectrum = np.split(np.copy(spectrum),2,axis=1)[1]
    
    (N,M) = spectrum.shape # Get number of points N and levels M
    if(cutoff): M = cutoff # Use cutoff if any
    part = np.zeros((N,M)) # Initialize partial current derivatives
    for i in range(N-1): # Loop over points
        for j in range(M): # Loop over levels
            part[i][j] += -(spectrum[i+1][j] - 2*spectrum[i][j] + spectrum[i-1][j]) / ((2*np.pi)/N)**2 # Calculate second-order derivative
    for j in range(M): # Update last point
        part[N-1][j] += -(spectrum[0][j] - 2*spectrum[N-1][j] + spectrum[N-2][j]) / ((2*np.pi)/N)**2
    total = np.sum(part, axis=1)
    return [np.transpose(part,axes=[1,0]), total]
    
#     [part, total] = supercurrent(spectrum, cutoff=cutoff)
#     return supercurrent_prime([part,total])

def supercurrent_prime(s=[None, None]): #
    [part, total] = s
    (M,N) = part.shape
    part_prime = [[] for m in range(M)]
    for i, p in enumerate(part):
        part_prime[i] = derivative_of_curve(p)
    total_prime = derivative_of_curve(total)
    return [part_prime, total_prime]

def critical_current_spec(spectrum, cutoff=None):
    [part, total] = supercurrent(spectrum, cutoff=cutoff)
    return critical_current([part, total])

def critical_current_from_current(s=[None, None]):
    [part, total] = s
    return np.amax(abs(total))

def sharpness_spec(spectrum, cutoff=None):
    [part, total] = supercurrent(spectrum, cutoff=cutoff)
    return sharpness_from_current([part, total])

def sharpness_from_current(s=[None, None]):
    [part, total] = s
    [part_prime,total_prime]=supercurrent_prime(s)
    s = np.log(abs(half_point(total_prime))/critical_current([part,total]))
    return s

def supercurrent_2d(spec):#spec.shape = (NPhi,NPhi,levels)
    N_phase = len(spec[:,:])
    curr_1, total_1 = [],[]
    for i in range(N_phase):
        (temp1,temp2) = supercurrent(spec[i])
        curr_1.append(temp1)
        total_1.append(temp2)
    curr_1 = np.asarray(curr_1)
    total_1 = np.asarray(total_1)

    curr_2, total_2 = [],[]
    for i in range(N_phase):
        (temp1,temp2) = supercurrent(spec[:,i])
        curr_2.append(temp1)
        total_2.append(temp2)
    curr_2 = np.transpose(np.asarray(curr_2),(1,0,2))
    total_2 = np.transpose(np.asarray(total_2),(1,0))
    return (curr_1, total_1, curr_2, total_2)

def supercurrent_3d(spec,b):#spec.shape = (NPhi,NPhi,field,levels)
    supercurrents=[[],[],[],[]]
    for field in range(len(b)):
        temp = supercurrent_2d(spec[:,:,field,:])
        for i in range(4): supercurrents[i].append(temp[i])
    for i in range(4): 
        supercurrents[i] = np.asarray(supercurrents[i])
    for i in (0,2):
        supercurrents[i] = np.transpose(supercurrents[i], (1,2,0,3))
    for i in (1,3):
        supercurrents[i] = np.transpose(supercurrents[i], (1,2,0))
    return supercurrents

