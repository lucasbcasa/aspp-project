"""
Collection of useful mathematical constants and common expressions.
"""
import numpy as np

#: Constant values that are often accessed by different modules
global hbar
global m
global ratio
global eff
global sigma_x
global sigma_y
global sigma_z
global sigmas

#: Planck's constant
hbar = 6.582119569*10**-13 # meV . s

#: Electron mass
m = ( 0.51099895*10E9 / ( (2.99792458*10E8)**2) ) * 10E-18 # meV . s^2 / nm^2

#: Ratio for the effective electron mass in InSb; Extracted from https://doi.org/10.1016/S1350-4495(99)00020-1
eff = 0.015 # unitless

#: Effective electron mass in InSb
m_eff = m*eff # meV . s^2 / nm^2

#: Useful factor that appears in the calculation of the hopping term
ratio_eff = hbar**2 / m_eff / 2 # meV . nm^2

#: Boltzmann constant
kB = 8.617333262145 * 10**(-5) # meV/mK

#: Array with superconducting phase values in radians
phi_space = np.linspace(0,2*np.pi,100,endpoint=False)

#: Array with superconducting phase values in units of pi
phi_space_pi = np.linspace(0,2,100,endpoint=False)

#: Array with k wavenumber values
k_space = np.linspace(-np.pi/2,np.pi/2,201)

#: Pauli matrices
sigma_x = np.array(((0, 1), (1, 0)))
sigma_y = np.array(((0, -1j), (1j, 0)))
sigma_z = np.array(((1, 0), (0, -1)))
sigmas = np.array((sigma_x, sigma_y , sigma_z))

#: Dictionary of parameters and their units
params_units = {'n_sc':' $[sites]$', 
                'n_normal':' $[sites]$', 
                'a':' $[nm]$',
                'delta':'$[meV]$',
                'alpha':'$[meV.nm]$',
                'mu':' $[meV]$',
                'mu_s': '$[meV]$',
                'mu_n': '$[meV$]',
                'B':' $[meV]$',
                'b':' $[B_c]$',
                'phi':' $[rad]$',
                'phi_pi':'$[\pi]$',
                'k':'$[a^{-1}]$',
                'kT':'$[meV]$',
                'T':'$[mK]$',
                'delta_mu':'$[meV]$'
               }

#: Dictionary of parameters and their plot labels
params_labels = {'n_sc':'$N_s~[sites]$', 
                 'n_normal':'$N_n~[sites]$', 
                 'a':'$a~[nm]$',
                 'delta':'$\Delta~[meV]$',
                 'alpha':'$alpha~[meV.nm]$',
                 'mu':'$\mu~[meV]$',
                 'mu_s':'$\mu_s~[meV]$',
                 'mu_n':'$\mu_n~[meV$]',
                 'B':'$B~[meV]$',
                 'b':'$B~[B_c]$',
                 'phi':'$\phi~[rad]$',
                 'phi_pi':'$\phi~[\pi]$',
                 'k':'$k~[a^{-1}]$',
                 'tau':'$tau$',
                 'kT':'$k_BT~[meV]$',
                 'T':'$T~[mK]$',
                 'delta_mu':'$\delta\mu~[meV]$'
                }

#: Dictionary of observables and their plot labels
observable_labels = {'absolute_sharpness': r'$\mathcal{S}_{abs}$',
                    'supercurrent_prime_pi': r'$-dI/d\phi_\pi$',
                    'spectra_normalized': r'$E/\Delta$',
                    'spectra': r'$E~[meV]$',
                    'supercurrent': r'$I~[meV~e/\hbar]$',
                     'supercurrent_prime': r'$-dI/d\phi~[meV~e/\hbar~rad]$',
                    'critical_current': r'$I_c~[meV~e/\hbar]$',
                    'critical_phi': r'$\phi_c~[rad]$'}

#: Alphabetical order to be used for labeling subplots
subplot_labels = list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
