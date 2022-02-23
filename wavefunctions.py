import numpy as np

def wavefunction(eig, N_probe=0):
    l = int(len(eig)/4)
    wave_nw = np.zeros(l)
    for i in range(l):
        wave_nw[i] = abs(eig[2*i])**2 + abs(eig[2*i+1])**2 + (+ abs(eig[2*(i+l)])**2 + abs(eig[2*(i+l)+1])**2)
    if (N_probe==0): return wave_nw
    return np.split(wave_nw, [l-N_probe])

def gamma(eig, *side): #Majorana transformation to obtain separate Majorana peaks; @Can be made more general
    sign = 1
    if(side): sign = side[0]
    l = int(len(eig)/4)
    wave = np.zeros(l)
    for i in range(l):
        wave[i] = abs(eig[2*i] + sign*eig[2*(i+l)])**2 + abs(eig[2*i+1] + sign*eig[2*(i+l)+1])**2
    return wave

def local(wave, r): # $Haven't used this in a long while; Maybe replace it with (exponential) curve-fitting function
    p=0
    for i in range(int(len(wave)*r)):
        p+=wave[i]
    return p