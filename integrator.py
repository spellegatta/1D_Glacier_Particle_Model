import numpy as np

from particle import Particle
from forces import gravity, collisions, normal_force, stiffness, buoyancy, fixed_constraint_force, constant_force, beam_force_1d



class Verlet_integrator:
    def __init__(self, dt):
        self.dt=dt


    def step(self, particles, close_particles, coll_coeff, gamma, k_suolo, young_mod, L0, A, bound_particles, constraints_bool, constant_force_bool, buoyancy_bool, w_level=None, w_density=None, glen_constant_force=None, time=None ):
        buoyancy_control= True if w_level==None or w_density==None else False
        Glen_control= True if constant_force==None or time==None else False
        if buoyancy_bool and buoyancy_control:
            print('Buoyancy is present, but density and/or level of the water was not specified')
            exit()
        if constant_force_bool and Glen_control:
            print('Constant force is present, but no value of force or time has been specified')
        for p in particles:
            a_old = p.f / p.mass
            p.coord += p.v*self.dt + 0.5*a_old*self.dt**2
            p.a_old = a_old   

        for p in particles:
            p.f_nuova = 0
            p.f_nuova+=gravity(p)
            p.f_nuova+=normal_force(particle=p, k_suolo=k_suolo)
            if buoyancy_bool:
                p.f_nuova+=buoyancy(w_level=w_level, w_density=w_density, particle=p)
                #print('galleggiamento',buoyancy(w_level=w_level, w_density=w_density, particle=p), 'gravità', gravity(p))
            #print('gravity: ', gravity(p), ' collisions: ', collisions(particles, close_particles, i, coll_coeff) , 'normal force: ', normal_force(particle=p, k_suolo=k_suolo))
        
        if constraints_bool:
            particles[0].f_nuova+=fixed_constraint_force(particle=particles[0], x0=particles[0].R, k=1E7)

        if constant_force_bool:
            constant_force(particle=particles[-1], c=glen_constant_force, time=time)
                
        
        collisions(particles, close_particles, bound_particles, coll_coeff, gamma=gamma) 
        stiffness(young_mod=young_mod, L0=L0, A=A, particles=particles, bound_particles=bound_particles) 

        for p in particles:
            a_new = p.f_nuova / p.mass
            p.v += 0.5*(p.a_old + a_new)*self.dt
            p.f = p.f_nuova

    def Glen_data_collection_forces(self, particles, beams, young_mod, A, L0):
        k = young_mod* A / L0

        expected_beams=[(i,i+1) for i in range(len(particles)-1)]
        forces=[]
        if beams!=expected_beams:
            raise Exception('At least one beam was broken during the simulation')
        for pi, pj in beams:
            fi, _ = beam_force_1d(particles[pi], particles[pj], k=k, L0=L0)
            forces.append(abs(fi))
        return forces
    
    def Glen_data_collection_velocities(self, particles, beams, young_mod, A, L0):
        expected_beams=[(i,i+1) for i in range(len(particles)-1)]
        dv_list=[]
        if beams!=expected_beams:
            raise Exception('At least one beam was broken during the simulation')
        else:
            for pi, pj in beams:
                dv_list.append(abs(particles[pi].v-particles[pj].v))

        return dv_list


                               


