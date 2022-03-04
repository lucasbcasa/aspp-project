import numpy as np
import time
import math
import scipy.optimize

import constants
import hamiltonians
import temperature
import sharpness

import warnings

# Finds the spectrum for a given set of parameters
def solve_single(nanowire, H=None, **kwargs):
    if H is None: 
        if kwargs.get('BdG', True): H = hamiltonians.HBdG(nanowire)
        else: H = hamiltonians.Hamiltonian(nanowire, k=nanowire.k)
    if kwargs.get('eigenvectors',False): 
        return np.linalg.eigh(H)
    else: 
        return np.linalg.eigvalsh(H)

# Finds the spectrum for a given set of parameters, varying only one of them (can't be a length!!!)
def solve(nanowire, quiet=True, **kwargs): # Uses only first set of parameters
    
    start = time.time()
    key = list(kwargs)[0]
    values = kwargs[key]
    
    steps = len(values)
    stops = np.linspace(0,9*steps//10,10, dtype=int)
    
    spectrum = np.empty((steps,4*nanowire.n_wire))
    for i, val in enumerate(values):
        setattr(nanowire, key, val)
        H = hamiltonians.HBdG(nanowire)
        spectrum[i] = np.linalg.eigvalsh(H)
        if not quiet and i in stops: print(round(100*i/steps,1),"%", end ='|')

    end = time.time()
    if not quiet: print("Time elapsed when solving for", nanowire.n_wire, "sites: ", round((end - start)/60,2), "min")
    return spectrum

# Finds the band spectrum for a given set of parameters
def solve_k(cell, k=constants.k_space, **kwargs):
    return solve(cell, k=k, **kwargs)

# Finds the spectrum for a given set of parameters, varying only field
def solve_field(junction, b=np.linspace(0,2,200,endpoint=False), **kwargs):
    return solve(junction, b=b, **kwargs)

# Finds the spectrum for a given set of parameters, varying only phase difference
def solve_phase(junction, phi=np.linspace(0,2*np.pi,200,endpoint=False), **kwargs):
    #print(junction.__dict__)
    return solve(junction, phi=phi, **kwargs)

# Takes in a junction object and a phase value and returns the total supercurrent at given phase
# Does not converge the value
# Adapted for temperature
def supercurrent_at(junction, phi, **kwargs):
    
    T = junction.T
    
    dphi = kwargs.get('dphi', 1E-4)
    junction.phi = phi + dphi
    spec_plus = np.split(np.linalg.eigvalsh(hamiltonians.HBdG(junction)), 2)[1]
    junction.phi = phi - dphi
    spec_minus = np.split(np.linalg.eigvalsh(hamiltonians.HBdG(junction)), 2)[1]
    
    junction.phi = phi
    spec_phi = np.split(np.linalg.eigvalsh(hamiltonians.HBdG(junction)), 2)[1] 
    
    # Spec_minus in front to counteract minus sign from equation
    supercurrent = temperature.smear_supercurrent(T,
                                                  spec_phi,
                                                  spec_minus - spec_plus) 
    
    return np.sum(supercurrent)/(2*dphi)

# Takes a junction object and returns the critical current
# Converges the value through bounded minimization procedure
# Not adapted for temperature
def critical_current(junction, **kwargs):
    min_opts = {'bounds':[0,np.pi], 'method':'bounded', 'options':{'xatol': 1E-4, 'maxiter': 500}}
    #min_opts = {'bounds':[[0,np.pi]], 'tol':1e-4}#, 'options':{'xatol': 1E-4, 'maxiter': 500}}
    res = scipy.optimize.minimize_scalar(lambda phi: -abs(supercurrent_at(junction, phi, **kwargs)), **min_opts)
    if kwargs.get('full_result',False): return res
    if (res.success): return abs(res.fun)
    else: return 'Maximization Failed'

# Takes in a junction object and a phase value and returns the total supercurrent derivative at given phase
# Does not converge the value
# Adapted for temperature
def supercurrent_prime_at(junction, phi, **kwargs):

    T = junction.T
        
    dphi = min(kwargs.get('dphi', 1E-4), abs(np.pi-kwargs.get('crit_phi',0))/2)
        
    if junction.phi!=phi: junction.phi = phi
    spec_phi = np.split(np.linalg.eigvalsh(hamiltonians.HBdG(junction)), 2)[1]
    junction.phi = phi + dphi
    spec_plus = np.split(np.linalg.eigvalsh(hamiltonians.HBdG(junction)), 2)[1]
    junction.phi = phi - dphi
    spec_minus = np.split(np.linalg.eigvalsh(hamiltonians.HBdG(junction)), 2)[1]
    
    supercurrent_prime = temperature.smear_supercurrent_prime(T, spec_phi, (spec_plus - spec_minus)/2, (spec_plus - 2*spec_phi + spec_minus))
    
    return np.sum(supercurrent_prime)/(dphi**2)
    

# Receives SNS_Junction object and returns the svalue of absolute sharpness;
# Not adapted for temperature
def absolute_sharpness_at(junction, phi, **kwargs):
    warnings.warn("Warning: Calling this method through the solvers library is deprecated. Call it from the sharpness library instead.")
    return sharpness.absolute_sharpness_at(junction, phi, **kwargs)

# Receives SNS_Junction object and returns the svalue of relative sharpness;
# Not adapted for temperature
def relative_sharpness_at(junction, phi, **kwargs):
    warnings.warn("Warning: Calling this method through the solvers library is deprecated. Call it from the sharpness library instead.")
    return sharpness.relative_sharpness_at(junction, phi, **kwargs)

# Receives SNS_Junction object and returns the value of absolute sharpness;
# Not adapted for temperature
def absolute_sharpness(junction, **kwargs):
    warnings.warn("Warning: Calling this method through the solvers library is deprecated. Call it from the sharpness library instead.")
    return sharpness.absolute_sharpness(junction, **kwargs)

# Receives SNS_Junction object and returns the value of relative sharpness;
# Not adapted for temperature
def relative_sharpness(junction, **kwargs):
    warnings.warn("Warning: Calling this method through the solvers library is deprecated. Call it from the sharpness library instead.")
    return sharpness.relative_sharpness(junction, **kwargs)

def test_func():
    print(2)
    return