# MASTER - MAcroscopic Simulation of Traffic to Enable Road predictions
----------------------------------------------------------------------

MASTER == **MA**croscopic **S**imulation of **T**raffic to **E**nable **R**oad predictions


## Description
--------------

MASTER is a non-GUI second-order macroscopic c/c++ traffic flow model
simulator, mainly for the Gas-Kinetic-based Traffic (GKT) model but
there are also other variants. Even diffusion-transport processes
(Fokker-Planck equations, Black Scholes option pricing can be
calculated. The simulator can also be "misused" to simulate
first-order (LWR-type) models.

Features:

- academic purpose for simulating the model dynamics on a single road
 with onramps, offramps, and bottlenecks. This is not a network
 traffic simulator

- multiple models though the main and most stable model is the GKT
  model. They  are selected by the choice_model entry of the .inp file
  of the corresponding project 

  * Gas-Kinetic-based Traffic (GKT) model in its original form:
    choice_model=0 

  * GKT model with a modification to make it more stable for high
    densities (choice_model=1)

  * Kuehne-Kerner-Konhaeuser (KKK) model (choice_model=2)

  * another secend-order model similar to the KKK model, the VMM model
    (choice_model=21) 

  * Black-Scholes model (choice_model=5)

  * macroscopic version of the Optimal Velocity Model with a
    triangular fundamental diagram 
    (choice_model=8, Tmin=Tmax and beta=0.0001, se ebelow)

  * macroscopic version of the Full Velocity Difference Model Optimal
    Velocity Model with a triangular fundamental diagram 
    (choice_model=8, Tmin=Tmax and beta>0, se ebelow)

 *  1d-Fokker-Planck equation (FPE) (choice_model=50)

 *  a variant of the FPE  (choice_model=51)




## Documentation
----------------

A mathematical description of the models as well as the basic concepts
can be found in the book [Traffic Flow
Dynamics](http://www.traffic-flow-dynamics.org) by Treiber/Kesting,
particularly the Chapter 9 (second-order macroscopic models, and also
Chapters 7 and 8 (general concepts and first-order macro models).

Documentation by example can be found in the [_sim_] directory 


## Installation
---------------

just type "make master"
which compiles the simulator and places an executable in ~/bin/


## Usage
--------

The simulations are organized in projects with several input and
output files (everything is text-based.
To run a simulation, type 

```bash
master projectName
```
in a command-line window. 

For example,
```bash
master ramps/ramps
```

will run the project *ramps* in the *ramps/* directory. Notice that
all output goes to the simulation directory, not the directory
*master* is called (best to call *master* in the corresponding
simulation directory)



## Demos
--------

Demo projects are in the *sim/* directory, ordered by topic.

Each project has a *.run* file which runs *master project*, and calls
the *gnuplot* program to make some graphical representation of the
results. If you do not want plotting, just call *master project*


## References 
-------------

[1] M. Treiber and A. Kesting. [Traffic Flow Dynamics, Data, Models and Simulation](http://www.traffic-flow-dynamics.org). Springer 2013. [Link](http://www.springer.com/physics/complexity/book/978-3-642-32459-8)

[2] D. Helbing, A. Hennecke, V. Shvetsov, and M. Treiber (2001)
MASTER: Macroscopic traffic simulation based on a gas-kinetic,
non-local traffic model. Transportation Research B 35, 183-211.

[3] M. Treiber, A. Hennecke, and D. Helbing (1999)
Derivation, Properties, and Simulation of a Gas-Kinetic-Based, Non-Local Traffic Model,
Phys. Rev. E 59, 239-253. 

[4] D. Helbing, A. Hennecke, and M. Treiber (1999)
Phase Diagram of Traffic States in the Presence of Inhomogeneities,
Phys. Rev. Lett. 82, 4360-4363. 

[5] M. Treiber and D. Helbing (1999)
Macroscopic Simulation of Widely Scattered Synchronized Traffic States
J. Phys. A: Math. Gen. 32, L17-L23.

[6] D. Helbing and M. Treiber (1999)
Numerical Simulation of Macroscopic Traffic Equations Computing in Science and Engineering (CiSE) 5, 89 (1999).

[7] D. Helbing and M. Treiber (1998)
Jams, waves, and clusters.,
Science 282, 2001-2003. DOI 10.1126/science.282.5396.2001 


