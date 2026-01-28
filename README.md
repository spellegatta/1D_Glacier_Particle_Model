# 1D Particle model for the simulation of glacier dynamics
In the notebook *Data_analysis.ipynb* the data calculated with the code in **run.py**, forces.py, integrator.py, particle.py and simulation.py will be displayed and the results will be shown. 

In the file particle.py it was implemented a class named *Particle*, which embeds the properties of each particle. 

In the file **forces.py** the forces acting on the particles were implemented. 

In the file **integrator.py** a Velocity Verlet integrator has been implemented in order to solve the equation of motion. 

In the file **simulation.py** the whole system of $N$ particles has been initialized.

In the **run.py** file, all these class were used in order to compute the simulation. In that code it is possible to change the parameters of the simulation and by running the command *python run.py* it is possible to make a new simulation.
