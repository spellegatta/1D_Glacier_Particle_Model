# 1D Particle model for the simulation of glacier dynamics
In the notebook *Data_analysis.ipynb* the data calculated with the code in **run.py**, forces.py, integrator.py, particle.py and simulation.py will be displayed and the results will be shown. 

In the file particle.py it was implemented a class named *Particle*, which embeds the properties of each particle. These properties are:
+ Its mass (each particle has a mass density of $900\pm 10 \, kg\cdot m^{-3}$)
+ Its side and radius (the particles have been treated throughout all the code as circles, but they have been initialized as hexagons)
+ Its coordinate
+ Its velocity
+ The forces in the past and present step (necessary for the Velocity Verlet integration)



In the file **forces.py** the forces acting on the particles were implemented. In the paper, the equation of motion was 
\begin{equation}
    \mathbf{M\ddot{r}}_i + \mathbf{C\dot{r}}_i+\sum_j\gamma_{ij} \mathbf{C}'\mathbf{\dot{r}}_{ij}+\sum_j \gamma_{ij}' \mathbf{Kr}_{ij}=F_i
\end{equation}
where $\mathbf{M}$ is the diagonal mass-matrix containing the masses and the moments of inertia, $\mathbf{r}_i$ and $\mathbf{\dot{r}}_i$ are the position and velocity of the particle $i$, $\mathbf{r}_{ij}$ and $\mathbf{\dot{r}}_{ij}$ are the corresponding relative vectors for particles $i$ and $j$. $\mathbf{C}$ is the matrix containing the damping coefficients for drag, while $\mathbf{C}'$ contains the coefficients for the inelastic collisions. The parameter $\gamma_{ij}$ is 0 for particles not in contact and 1 ptherwise, while $\gamma_{ij}'$ is unity ofr connected particles and 0 otherwise. $\mathbf{K}$ is the stiffness matrix and $F_i$ is the sum of other forces acting on the particle.

In this code a slightly different approach than the one used in the paper was used. This was done because for a 1-dimensional simulation, using the same approach that was used for the 2-dimensional problem meant overcomplicating the problem. The forces that acted on the particle in the 1D simulation are:
+ gravity: $-\rho_i Vg$.
+ collisions: if the particles $i$ and $j$ collided the force acting on the $i$-th particle was $F_i=\gamma\mathbf{\Delta x}-C'\Delta v$, where $C'=10^5 \, Nsm^{-1}$, $\mathbf{\Delta x}$ is the overlap of the 2 particles and is written as a vector because it contains the sign that makes the force act in the opposite direction of the collision, $\gamma = 10^6\, Nm^{-1}$ (chosen empirically). $\Delta v$ is the relative velocity. The force acting on the $j$-th particle was $-F_i$ according to Newton's third law. The term depending on the overlap has been added in this model, in order not to make one particle penetrate another. It has not been made clear but in the paper they probably used a constraint to do so.
+ normal force: it has not been implemented as a constraint, but as an elastic force added to a velocity damping. This was done to simulate a muddy floor as was done in one of the simulations of the paper. 
+ Interaction potential between two particles: in the paper it was modelled as en Euler-Bernoulli beam. In 1 dimension that meant that angular momenta and inertia should have been included, complicating by a great deal the code. Therefore, another model suggested by the paper was chosen. It was used an harmonic potential between the particles with respect to their positions, neglecting the harmonic potential with respect to their node rotations away from the axis, which would have seemed rather useless in a 1D model.

In the file **integrator.py** a Velocity Verlet integrator has been implemented in order to solve the equation of motion. This integrator update the positions using current velocities and accelerations, reset the forces and computes them at the new positions and updates velocities using the average of old and new accelerations.

In the file **simulation.py** the whole system of $N$ particles has been initialized. Other than store the particles and their characteristics, the class *system* store the bonds between the particles and the bonds that were broken. Other than the function that break the beams that accumulated too much elastic energy (the elliptical criterion was not implemented here, due to its impracticality in a 1D model) and the function that computes which particles are close for the computation of collision forces, the class contains function as *stochastic_melting* and *stochastic_refreezing* needed to be able to simulate viscous behaviour.

In the **run.py** file, all these class were used in order to compute the simulation. In that code it is possible to change the parameters of the simulation and by running the command *python run.py* it is possible to make a new simulation.
