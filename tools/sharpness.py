import solvers
import numpy as np

# Receives SNS_Junction object and returns the svalue of absolute sharpness;
def absolute_sharpness_at(junction, phi, **kwargs):
    curr_prime = kwargs.get('curr_prime', solvers.supercurrent_prime_at(junction, phi, **kwargs)) # Get or calculate the current prime at phi
    if curr_prime==0: return math.inf
    return np.log( abs( curr_prime ) )

# Receives SNS_Junction object and returns the svalue of relative sharpness;
def relative_sharpness_at(junction, phi, **kwargs):
    crit_curr = kwargs.get('crit_curr', solvers.critical_current(junction,**kwargs)) # Get or calculate the critical current
    curr_prime = kwargs.get('curr_prime', sovers.supercurrent_prime_at(junction, phi, **kwargs)) # Get or calculate the current prime at phi
    if curr_prime==0: return math.inf
    return np.log( abs( curr_prime / crit_curr ) )

# Receives SNS_Junction object and returns the value of absolute sharpness;
def absolute_sharpness(junction, **kwargs):
    if 'curr_prime' in kwargs:
        return np.log( abs( kwargs['curr_prime'] ) )
    return solvers.absolute_sharpness_at(junction, np.pi, **kwargs)

# Receives SNS_Junction object and returns the value of relative sharpness;
def relative_sharpness(junction, **kwargs):
    if 'curr_prime' in kwargs and 'crit_curr' in kwargs:
        return np.log( abs( kwargs['curr_prime'] / kwargs['crit_curr'] ) )
    return solvers.relative_sharpness_at(junction, np.pi, **kwargs)