import numpy as np
import constants

def inhomogenousPotential(nanowire, delta_mu, omega = 0.05):
    if not hasattr(nanowire, 'n_sc'):
        n = nanowire.n_wire
        mu_s = nanowire.mu
        mu = mu_s - delta_mu * (np.exp(-0.5 * (np.array(range(n))**2) / (omega*n)**2 ) +
                         np.exp(-0.5 * ((n-1) - np.array(range(n)))**2 / (omega*n)**2 ) )
        return mu
    else:
        n_sc = nanowire.n_sc
        n_normal = nanowire.n_normal
        mu_s = nanowire.mu_s
        mu_s_profile = mu_s - delta_mu * (np.exp(-0.5 * (np.array(range(n_sc))**2) / (omega*n_sc)**2 ) + 
                                          np.exp(-0.5 * ((n_sc-1) - np.array(range(n_sc)))**2 / (omega*n_sc)**2 ) )
        mu = np.zeros(nanowire.n_wire)
        mu[range(n_sc)] = mu_s_profile # Left S region
        mu[range(n_sc, n_sc + n_normal)] = np.full((n_normal), nanowire.mu_n) # N region
        mu[range(n_sc + n_normal, 2*n_sc + n_normal)] = mu_s_profile # Right S region
        
        return mu

class Nanowire:
    def __init__(self, **kwargs):
        self.n_wire = kwargs.get('n_wire', 100) #('key', default value) -> creates the entry as default value if no value is found
        
        self._a = kwargs.get('a', 10) # nm
        
        self._t = constants.ratio_eff / (self._a**2) # meV
        
        self.alpha = kwargs.get('alpha', 20) # nm.meV
        self._tso = self.alpha/(2*self._a)
        self.so_axis = kwargs.get('so_axis', (0,1,0)) # Spin-orbit axis
        
        self._mu = kwargs.get('mu', 0.5) # meV
        
        self._delta = kwargs.get('delta', 0.25) # meV
        
        self.Bc = kwargs.get('Bc', np.sqrt(self._delta**2 + self._mu**2)) # meV
        
        self.b = kwargs.get('b', 0) # units of Bc
        self.B = kwargs.get('B', self.b*self.Bc) #meV
        
        self.periodic_bc = kwargs.get('periodic_bc', False)
        self.k = kwargs.get('k', None)
        
        self.T = kwargs.get('T', 0)
        
    @property
    def a(self): return self._a
    @a.setter
    def a(self, a):
        self._a, self._t = a, constants.ratio_eff / (a**2)
        self._tso = self.alpha/(2*self._a)
        
    @property
    def t(self): return self._t
    @t.setter
    def t(self, t):
        self._t, self._a = t, np.sqrt(constants.ratio_eff / t)
        self._tso = self.alpha/(2*self._a)
        
    @property
    def tso(self): return self._tso
    @tso.setter
    def tso(self, tso):
        self._tso, self._a = tso, self.alpha/(2*tso)
    
    
    @property 
    def delta(self): return self._delta
    @delta.setter 
    def delta(self, delta): 
        self._delta, self.Bc = delta, np.sqrt(delta**2 + self.mu**2)
        
    @property
    def delta_profile(self):
        try: return self._delta_profile
        except AttributeError: self.build_profile( _delta_profile=[[self.n_wire, self.delta]] )
        return self._delta_profile
    @delta_profile.setter
    def delta_profile(self, delta_profile):
        self._delta_profile  = delta_profile

# Chemical Potential        

    @property 
    def mu(self): return self._mu
    @mu.setter 
    def mu(self, mu): 
        self._mu, self.Bc = mu, np.sqrt(self.delta**2 + mu**2)
        self.build_profile( _mu_profile=[[self.n_wire, mu]] )
        
    @property
    def mu_profile(self):
        try: return self._mu_profile
        except AttributeError: self.build_profile( _mu_profile=[[self.n_wire, self.mu]] )
        return self._mu_profile
    @mu_profile.setter
    def mu_profile(self, mu_profile):
        self._mu_profile = mu_profile
        
    #  Inhomogenous chemical potential: if delta_mu is set, then create profile   
    @property
    def delta_mu(self): return self._delta_mu
    @delta_mu.setter
    def delta_mu(self, delta_mu):
        self._delta_mu = delta_mu
        self._mu_profile = inhomogenousPotential(self, delta_mu)

# Zeeman Field

    @property 
    def B(self): return self._B
    @B.setter 
    def B(self, B): self._B, self._b = B, B / self.Bc
        
    @property 
    def b(self): return self._b
    @b.setter 
    def b(self, b): self._b, self._B = b, b * self.Bc
        
    #@ need to implement constinuous transitions
    def build_profile(self, **kwargs): # Takes key=pieces;
                                       # piece = [length, value] 
        for key in kwargs:
            attribute = np.array([])
            for piece in kwargs[key]:
                attribute = np.concatenate((attribute, np.full((piece[0]), piece[1])))
            self.__dict__[key] = attribute
    
    def update_profile(self, **kwargs): # Takes key=pieces;
                                        # piece = [x0, length, new_val]
        for key in kwargs:
            attribute = self.__dict__[key]
            for piece in kwargs[key]:
                attribute[piece[0]:piece[0]+piece[1]] = piece[2]

class SNS_Junction(Nanowire):
    def __init__(self, **kwargs):
        Nanowire.__init__(self, **kwargs)

        # Set attributes directly because n_wire can't be calculated yet
        self._n_left = kwargs.get('n_left', kwargs.get('n_sc', 100))
        self._n_right = kwargs.get('n_right', kwargs.get('n_sc', self._n_left))
        self._n_normal = kwargs.get('n_normal', 4)
        
        self.update_n_wire()

        self._phi = kwargs.get('phi', 0)
        self.build_delta_profile()
        
        self._mu_s = kwargs.get('mu_s', self._mu)
        self._mu_n = kwargs.get('mu_n', self._mu)
        self.build_mu_profile()
        
        self.tau = kwargs.get('tau', 1)
    
    def build_delta_profile(self):
        pieces = [[self.n_left, self.delta * np.exp(1j * self.phi)], [self.n_normal, 0], [self.n_right, self.delta]]
        self.build_profile( _delta_profile=pieces )
    
    # Updates always the left SC
    def update_delta_profile(self):
        self.update_profile( _delta_profile=[[0, self.n_left, self.delta * np.exp(1j * self.phi)]] )
        
    def build_mu_profile(self):
        pieces = [[self.n_left, self.mu_s], [self.n_normal, self.mu_n], [self.n_right, self.mu_s]]
        self.build_profile( _mu_profile=pieces )
        
    def update_mu_profile_s(self):
        self.update_profile( _mu_profile=[[0, self.n_left, self.mu_s], [self.n_left + self.n_normal, self.n_right, self.mu_s]] )
        
    def update_mu_profile_n(self):
        self.update_profile( _mu_profile=[[self.n_left, self.n_normal, self.mu_n]] )
        
    def update_n_wire(self):
        self.n_wire = self._n_left + self._n_normal + self._n_right
    
    @property
    def n_sc(self): 
        if self._n_left==self._n_right: return self._n_left
        else: return (self._n_left, self._n_right)
    @n_sc.setter
    def n_sc(self, n_sc):
        self._n_left = n_sc
        self._n_right = n_sc
        self.update_n_wire()
        self.build_delta_profile()
        self.build_mu_profile()
    
    @property
    def n_left(self): return self._n_left
    @n_left.setter
    def n_left(self, n_left):
        self._n_left = n_left
        self.update_n_wire()
        self.build_delta_profile()
        self.build_mu_profile()

    @property
    def n_right(self): return self._n_right
    @n_right.setter
    def n_right(self, n_right):
        self._n_right = n_right
        self.update_n_wire()
        self.build_delta_profile()
        self.build_mu_profile()

    @property
    def n_normal(self): return self._n_normal
    @n_normal.setter
    def n_normal(self, n_normal):
        self._n_normal = n_normal
        self.update_n_wire()
        self.build_delta_profile()
        self.build_mu_profile()

    @property
    def phi(self): return self._phi
    @phi.setter
    def phi(self, phi):
        self._phi = phi
        self.update_delta_profile()
        
    # Phase in units of Pi
    @property
    def phi_pi(self): return self._phi / np.pi
    @phi_pi.setter
    def phi_pi(self, phi_pi):
        self._phi = phi_pi * np.pi
        self.update_delta_profile()
        
    @property
    def mu_s(self): return self._mu_s
    @mu_s.setter
    def mu_s(self, mu_s):
        self._mu_s = mu_s
        self.update_mu_profile_s()
        
    @property
    def mu_n(self): return self._mu_n
    @mu_n.setter
    def mu_n(self, mu_n):
        self._mu_n = mu_n
        self.update_mu_profile_n()
        
    #@property
    #def mu(self): return self.mu_profile
    #@mu.setter
    #def mu(self, mu):
    #    self._mu_n = mu
    #    self._mu_s = mu
    #    self.build_mu_profile()

def main():
    nanowire = Nanowire()

if __name__ == "__main__":
    main()
