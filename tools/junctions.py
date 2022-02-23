import numpy as np
import constants

import nanowires

class Junction(Nanowire):
    def __init__(self, **kwargs):
        Nanowire.__init__(self, **kwargs)
        
        self._geometry = kwargs.get('geometry') #[n_sc1, n_normal1, n_sc2, n_normal2,...]
        
        #Build delta profile from geometry

def main():
    nanowire = Nanowire()

if __name__ == "__main__":
    main()
