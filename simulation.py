import numpy as np
import math
import pandas as pd
import os 

from particle import Particle
from forces import gravity, collisions, normal_force
from integrator import Verlet_integrator

class system:
    def __init__(self, nparticles, side_lenght, typical_density, density_std):
        self.side_lenght=side_lenght
        self.R=self.side_lenght*pow(3, 0.5)/2
        self.typical_density=typical_density
        self.density_std=density_std
        self.n=nparticles
        self.particles=[]
        self.beams=[(i,i+1) for i in range(self.n-1)]
        self.broken_beams=[]

    def build(self):
        for j in range(self.n):
            x=2*(j)*self.R+self.R
            p=Particle(side_lenght=self.side_lenght, typical_density= self.typical_density, density_std=self.density_std, x=x)
            self.particles.append(p)

    def build_fall(self):
        for j in range(self.n):
            x=2.2*j*self.R+self.R
            p=Particle(side_lenght=self.side_lenght, typical_density= self.typical_density, density_std=self.density_std, x=x)
            self.particles.append(p)

    def out_positions(self, path, time):
        x=[]
        time_list=[time]*len(self.particles)
        for particle in self.particles:
            x.append(particle.coord)
        
        file_nuovo=not os.path.exists(path)
    
        df=pd.DataFrame({"time":time_list, "x":x})
        df.to_csv(path_or_buf=path, index=False, mode='a', header=file_nuovo)

    def close_particles(self, threshold):
        close_particles=[]
        for i, pi in enumerate(self.particles):
            for j, pj in enumerate(self.particles):
                if i>j:
                    d=abs(pi.coord-pj.coord)
                    if d<2*threshold*pi.R:
                        close_particles.append((i,j))
        return close_particles
    
    def fractures(self, critical_dx, young_mod, A, L0):
        k=young_mod*A/L0
        U_crit=0.5*k*critical_dx**2
        to_break=[]
        for (pi, pj) in self.beams:
            #print('pi, pj', pi, pj)
            U=0.5*k*(self.particles[pj].coord-self.particles[pi].coord-L0)**2
            if U>U_crit:
                to_break.append((pi, pj))
        for beam in to_break:
            self.beams.remove(beam)
            self.broken_beams.append(beam)
    
    def stochastic_melting(self, critical_dx, young_mod, A, L0, U_frac, hom_T, m_arr, dt):
        k=young_mod*A/L0
        U_crit=0.5*k*critical_dx**2
        to_break=[]
        Arrhenius = m_arr * 2.95 * 1E-9 * math.exp(-7.88*1E4 / (8.321 * hom_T) + 3 * 0.17 / (273.39 - hom_T)**1.17)
        for (pi, pj) in self.beams:
            #print('pi, pj', pi, pj)
            U=0.5*k*(self.particles[pj].coord-self.particles[pi].coord-L0)**2
            lamda = 2 * Arrhenius * U * young_mod ** 2 /(A**2)
            P=1-math.exp(-lamda*dt)

            if U>(U_frac*U_crit) and np.random.rand()<P:
                to_break.append((pi, pj))
        for beam in to_break:
            self.beams.remove(beam)
            self.broken_beams.append(beam)


    def stochastic_refreezing(self, critical_dx, young_mod, A, L0, U_frac, hom_T, m_arr, dt, refreezing_fraction):
        k=young_mod*A/L0
        U_crit=0.5*k*critical_dx**2
        to_create=[]
        Arrhenius = m_arr * 2.95 * 1E-9 * math.exp(-7.88*1E4 / (8.321 * hom_T) + 3 * 0.17 / (273.39 - hom_T)**1.17)
        #print(Arrhenius)
        lamda = 2 * Arrhenius * U_crit * young_mod ** 2 /(A**2)
        mu=refreezing_fraction*lamda
        P=1-math.exp(-mu*dt)
        
        for (i, j) in self.broken_beams:
            #print('pi, pj', pi, pj)
            U=0.5*k*(self.particles[j].coord-self.particles[i].coord-L0)**2


            if U<(U_frac*U_crit) and np.random.rand()<P:
                to_create.append((i, j))
        for beam in to_create:
            self.beams.append(beam)
            self.broken_beams.remove(beam)

    def out_Glen(self, force, A, dv, L0, path):

        file_nuovo=not os.path.exists(path)
        df=pd.DataFrame({"Stress":force/A, "strain": dv/L0})
        df.to_csv(path_or_buf=path, index=False, mode='a', header=file_nuovo)







 
            

