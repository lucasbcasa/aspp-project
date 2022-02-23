import numpy as np
import constants

# Takes in temperature, an array of energies and an array of the supercurrent contributions in a fixed phase difference
def smear_supercurrent(T, energy, supercurrent, ):
    if T==0:
        return supercurrent # If the temperature is zero do nothing
    kT = constants.kB * T
    return np.tanh( energy / (2*kT) ) * supercurrent # If the temperature is finite smear the currents
    
# Takes in temperature, array of energies, array of supercurrent 
# and array of supercurrent_prime contributions in a fixed phase difference
def smear_supercurrent_prime(T, energy, supercurrent, supercurrent_prime):
    if T==0:
        return supercurrent_prime # If the temperature is zero do nothing
    # Taking the derivative of the smeared supercurrent with respect to the phase difference and using the product rule:
    kT = constants.kB * T
    return  ( (1 - np.tanh( energy / (2*kT) )**2) * ( supercurrent / (2*kT) ) * supercurrent ) + (np.tanh( energy / (2*kT) ) * supercurrent_prime)