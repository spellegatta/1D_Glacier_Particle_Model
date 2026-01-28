import numpy as np

class Particle:
    def __init__(self, side_lenght, typical_density, density_std, x):
        self.coord=x
        self.v=0.
        self.f=0.
        self.f_nuova=0.
        self.side=side_lenght
        self.R=side_lenght*np.pow(3, 0.5)/2
        self.mass=np.random.normal(typical_density, density_std)*np.pi*self.R**2


        