import os
from simulation import system
from integrator import Verlet_integrator
import numpy as np

sim=system(nparticles=10, side_lenght=1, typical_density=900, density_std=10)
sim.build()

dt=1e-4
n_steps_partial=10
critical_dx=2.1*sim.particles[0].R
E=1E5
A=np.pi*(sim.particles[0].R**2)
L0=(2*sim.particles[0].R)
U_frac = 1e-1
hom_T = 272
m_arr=8
refreezing_fraction=1e-3
path='data/positions_buoyancy_2.csv'
buoyancy_bool=True
w_level= 15**0.5 * 5
w_density=1025


integr=Verlet_integrator(dt=dt)
for i in range(5000000):
    integr.step(particles=sim.particles, close_particles=sim.close_particles(threshold=1), coll_coeff=1E6, gamma=1E5, k_suolo=1E5, young_mod=E, L0=L0, A=A, bound_particles=sim.beams, buoyancy_bool=buoyancy_bool,constant_force_bool=False, constraints_bool=False ,w_level=w_level, w_density=w_density)
    sim.fractures(critical_dx=critical_dx, young_mod=E, A=A, L0=L0)
    sim.stochastic_melting(critical_dx=critical_dx, young_mod=E, A=A, L0=L0, U_frac=U_frac, hom_T=hom_T, m_arr=m_arr, dt=dt)
    sim.stochastic_refreezing(critical_dx=critical_dx, young_mod=E, A=A, L0=L0, U_frac=U_frac, hom_T=hom_T, m_arr=m_arr, dt=dt, refreezing_fraction=refreezing_fraction)
    #print(sim.beams)
    if os.path.exists(path) and i==0:
        os.remove(path)
    
    if (i%n_steps_partial)==0:
        time=n_steps_partial*dt*i
        sim.out_positions(path=path, time=time)
