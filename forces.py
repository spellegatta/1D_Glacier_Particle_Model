import numpy as np
from scipy import constants
import math

from particle import Particle

def gravity(particle):
    gravity=particle.mass*constants.g
    return -gravity

def collisions(particles, close_particles, bound_particles, coll_coeff, gamma):

    for (p1, p2) in close_particles:
        if (p1, p2) in bound_particles:
            continue
        dv = particles[p1].v - particles[p2].v
        overlap = particles[p1].R + particles[p2].R - abs(particles[p1].coord - particles[p2].coord)
        if overlap <= 0:
            continue
        sign=1 if particles[p2].coord>particles[p1].coord else -1
        particles[p1].f_nuova-=coll_coeff*sign*overlap+gamma*dv
        particles[p2].f_nuova+=coll_coeff*sign*overlap+gamma*dv
        
        #print('distanze e forze', (overlap, collision_force))
            

def normal_force(particle, k_suolo):


    if particle.coord < particle.R:
        gamma = 2 * np.sqrt(particle.mass * k_suolo)

        overlap = particle.R - particle.coord

        F = k_suolo * overlap

        F -= gamma * particle.v

        #print('parte elastica: ', k_suolo * overlap, 'smorzamento: ', - gamma * particle.v)
        return F

    return 0.0


def beam_force_1d(pi, pj, k, L0):
    m_eff=pi.mass*pj.mass/(pi.mass+pj.mass)
    gamma = 2 * (m_eff * k)**0.5

    xi = pi.coord
    xj = pj.coord
    vi = pi.v
    vj = pj.v

    dx = xj - xi

    f = k * (dx - L0) - gamma*(vi-vj)
    return f, -f

def stiffness(young_mod, L0, A, particles, bound_particles):
    
    k = young_mod* A / L0

    for (i,j) in bound_particles:
        fi,fj=beam_force_1d(particles[i], particles[j], k, L0)
        particles[i].f_nuova+=fi
        particles[j].f_nuova+=fj
    

def buoyancy(w_level, w_density, particle):
    if (particle.coord+particle.R)<w_level:
        V=np.pi*particle.R**2
    elif (particle.coord-particle.R)>w_level:
        V=0
    else:
        if particle.coord<w_level:
            ratio = (w_level-particle.coord)/particle.R
            if ratio>1 and ratio<1.1:
                ratio=1
            theta= 2*math.acos(ratio)
            A_piccola= 0.5* particle.R**2 * (theta-math.sin(theta))
            V= np.pi*particle.R**2 - A_piccola
        else:
            ratio = (w_level-particle.coord)/particle.R

            if ratio<-1 and ratio>-1.1:
                ratio=-1
            theta= 2*math.acos(ratio)
            V= 0.5* particle.R**2 * (theta-math.sin(theta))
    return w_density*V*constants.g

def fixed_constraint_force(particle, x0, k):
    dx = particle.coord - x0
    gamma = 2 * np.sqrt(particle.mass * k)
    return -k * dx - gamma * particle.v

def constant_force(particle, c, time):
    particle.f_nuova+=c#*time




