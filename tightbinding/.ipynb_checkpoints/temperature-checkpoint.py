"""
This module processes the partial supercurrents data in order to implement temperature effects into it.

It can generate finite temperature data for both supercurrent and zero-frequency susceptibility profiles.
"""
import numpy as np
import constants

# Takes in temperature, an array of energies and an array of the supercurrent contributions in a fixed phase difference
def smear_supercurrent(T, energy, supercurrent):
    """
    Add temperature effects to the supercurrent contribution of a given state for a fixed phase value. Can process one or several energy levels at once.
    
    Parameters
    ----------
    T : float
        Temperature of the system, in Kelvin. It is then converted to meV through the Boltzmann constant.
    
    energy : float, 1-d np.array
        Energy for the corresponding state(s) at the fixed phase. It is necessary because the supercurrent contributions get suppressed depending on their energy to temperature ratio.
        
    supercurrent : float, 1-d np.array
        The state's contribution to the supercurrent at zero temperature. Array length must match the `energy` one.
        
    Returns
    -------
    float, 1-d np.array
        The state's contribution to the supercurrent, after temperature effect is applied. If 
    """
    if T==0:
        return supercurrent # If the temperature is zero do nothing
    kT = constants.kB * T # Get the temperature in units of energy
    return np.tanh( energy / (2*kT) ) * supercurrent # If the temperature is finite smear the currents
    
# Takes in temperature, array of energies, array of supercurrent 
# and array of supercurrent_prime contributions in a fixed phase difference
def smear_supercurrent_prime(T, energy, supercurrent, supercurrent_prime):
    """
    Add temperature effects to the contribution to the zero-frequency susceptibility of a given state for a fixed phase value.
    
    Parameters
    ----------
    T : float
        Temperature of the system, in Kelvin. It is then converted to meV through the Boltzmann constant.
    
    energy : float, 1-d np.array
        Energy for the corresponding state(s) at the fixed phase. It is necessary because the supercurrent contributions get suppressed depending on their energy to temperature ratio, and so the effects on the zero-frequency susceptibility also depend on this ratio.
        
    supercurrent : float, 1-d np.array
        The state's partial contribution to the supercurrent at zero temperature. Array length must match the `energy` one.
        
    supercurrent_prime : float, 1-d np.array
        The phase derivative of the state's contribution to the supercurrent at zero temperature. Array length must match the `energy` one.
        
    Returns
    -------
    float, 1-d np.array
        The state's contribution to the zero-frequency susceptibility, after temperature effect is applied.
    """
    if T==0:
        return supercurrent_prime # If the temperature is zero do nothing
    # Taking the derivative of the smeared supercurrent with respect to the phase difference and using the product rule:
    kT = constants.kB * T # Get the temperature in units of energy
    return  ( (1 - np.tanh( energy / (2*kT) )**2) * ( supercurrent / (2*kT) ) * supercurrent ) + (np.tanh( energy / (2*kT) ) * supercurrent_prime)