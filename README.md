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
  of the corresponding project. Examples are the following.

  * Gas-Kinetic-based Traffic (GKT) model in its original form (`choice_model=0`)  with some modifications to make it more stable for high densities. See [Traffic Flow Dynamics](http://www.traffic-flow-dynamics.org) by Treiber/Kesting, second edition, 2026.
     
  * Kuehne-Kerner-Konhaeuser (KKK) model (`choice_model=2`)

  * Black-Scholes model (`choice_model=5`)

  * macroscopic version of the Optimal Velocity Model and the FVDM with different 
    fundamental diagrams and an optional "three-phase" indifference zone (`choice_model=8`) 


  *  1d-Fokker-Planck equation (FPE) (choice_model=50) and a variant thereof

 




## Documentation
----------------

A mathematical description of the models as well as the basic concepts
can be found in the book [Traffic Flow
Dynamics](http://www.traffic-flow-dynamics.org) by Treiber/Kesting, second edition, Springer, 2026
particularly the Chapter 10 (second-order macroscopic models, and also
Chapters 8 and 9 (general concepts and first-order macro models).

Documentation of the actual simulation projects can be found in the [_Demos_] section.


## Installation
---------------

just type 

```
make master
```

which compiles the simulator and places an executable in *~/bin/*


## Running the program
----------------------

The simulations are organized in projects with several input and
output files (everything is text-based.
To run a simulation, type 

```bash
master projectName
```
in a command-line window. 

For example,
```bash
master sim/ramps
```

will run the project *sim/ramps* in the *sim/* directory. Notice that
all output goes to the simulation directory, not the directory
*master* is called from (best to call *master* in the corresponding
simulation directory, though)



## Demos
--------

Demo projects are in the *sim/* directory, ordered by topic.

Each project has a *.run* file which runs *master project*, and calls
the *gnuplot* program to make some graphical representation of the
results. If you do not want plotting, just call *master project*.
Each project has mandatory and optional files. Comment lines begin with a '%'.

Let's assume that the project is called `proj`

## Mandatory project file: `proj.inp`
-------------------------------------

Main control file. The example files are essentially self-explaining. Points to watch out for:

* Lines beginning with `%` are comments

* The model parameters (first 7 data lines) are for the GKT and its variants. They are ignored and/or overridden by optional parameters in other (conditionally mandatory) files if other macromodels are simulated. Also the GKT parameters are overridden as a function of following optional files:

  - `proj.Tr` or `proj.v0` override v0 and T by space-dependent values, respectively, allowing to model flow-conservative bottlenecks (see below),

  - `proj.A`: The assumed relative speed variance A(rho) (squared variation coefficient) of the GKT model is a tanh-like function of the density rho. Its amplitude parameters A0=A(0) and dA=(A(rhomax)-A0)/2 are in the `.inp` file while its positional parameters pos_rel=0.27 (relative density of maximum increase) and drho=0.1 (relative width of the increase) have default values in the code but can be overridden by the optional .A file   

* The following 7 data lines are intitial conditions. If not overridden by optional files and `choice_init`>2 (below), they denote the initial density and a possible localized initial perturbation or shock and their specifications like amplitude and position of the perturbation. For `choice_init`=1 or 2, the speed is initially equal to the equilibrium speed.

* The line with `choice_IC` selects the initial conditions. The values are

  - 0: periodic: rho = rhoinit + ampl_init*periodic(wavel_init) * Gauss envelope(w_init), Q=Qe(rho)
  - 1: shocks:  rho = rhoinit +(x in [x1,x2]) ? ampl_init : 0, x1/2=xmax*center_init +/- w_init/2, Q=Qe(rho)
  - 2: rho: Asymmetric localized dipole-like perturbation a la Kerner, Q=Qe(rho)
  - 3: rho from `proj.IC`, Q=Qe(rho)
  - 4: rho and Q from `proj.IC`

* the line with `choice_BC` selects the boundary conditions, possibly using the files `proj.BCl` and `proj.BCr`. The values are

  - 0: periodic BC
  - 1: Upstream fixed density, downstream free BC
  - 2: Upstream fixed density, downstream homogeneous von Neumann BC
  - 3: Upstream and downstream fixed (Dirichlet) densities, steady-state flows
  - 4: Upstream and downstream free BC
  - 5: Upstream and downstream prescribed densities and flows
  - 6: Upstream and downstream homogeneous von Neumann BC
  - 7: Upstream density an flow, downstream von Neumann
  - 8: After an intelligent determination whether traffic flow is free or congested at the boundaries, flow and densiy is used for free upstream and congested downstream situation. Otherwise, homogeneous von Neumann BCs are used.
  - 9: If upstream is free traffic or downstream congested, only the FLOWS are used while the densities are the inverse of the respective branch of the steady-state flow-density relations (very effective!), otherwise,  homogeneous von Neumann BCs are used.
  - 10: Upstream as for choice_BC==8), downstream, if congested, the speed (or Q/rho) is used while flow and density are calculated by (the correct branch of) the steady-state relations. If free, homogeneous downstream von Neumann BCs are used.

The different specifications mean the following:

  - Homogeneous von Neumann: zero gradient at the boundary
  - free: unchanged gradient at the boundary (reflecting minimal BC influence)

* The line setting `choice_rmp` specifies the network as follows: `choice-rmp`= ...

  - 0: No ramps
  - 1: A single ramp specified by `proj.rmp` (obsolete, overseeded by `choice_rmp`=2)
  - 2: One or several ramps whose center merging locations and merging lengthd are specified in `proj.rmps` while the flow of the individual ramps are specified in the files `proj.rmp<n>` with n=1,2,... (a negative flow denotes an off-ramp)
  - 3: An on-ramp with controlled inflow
  - 4: Virtual ramps denoting changes in the number of lanes by relative changes in the flow per lane (the whole MASTER simulator reduces multiple lanes to an effective single lane). For example
    * from 1 to 2 lanes: Qrmp=-0.5  (50% reduction per lane)
    * from 2 to 1 lane: Qrmp=1 (100% increase per lane)
    * From 3 to 2 lanes: Qrmp=0.5 (50% increase per lane)


* Block for numerical parameters

Self-explaining. It is up to the user to check the CFL conditions such as dx/dt<maximum propagation speed or D*dt/dx^2<1/2. Besides the unavoidable numerical diffusion in first-order schemes, a small explicit diffusion can be used to stabilize the model with respect to relaxational instabilities for very high densities (very small conditional diffusions of a maximum of 10 m^2/s are built into the code, even for D1=D2=0 in the `.inp` file to ensure stability) 

* Block for the output

 - dtout and dxout defines the output grid for the spatiotemporal density and flow written in `proj.dat`
 
 - if `choice_outp`=1, additional time and space slices are written to `proj.x<value>`, `proj.t<value>` as prescribed by the optional files `proj.t` and `proj.x`, respectively

* Choice of the macroscopic model: `choice_model`== ...

  - 0: GKT model (best use the upwind scheme, `choice_method`=1)
  - 1: GKT model with high density correction (obsolete, now built into the standard model)
  - 2: Kerner-Konhaeuser model (requires `choice_method`=0)
  - 3: Macroscopic version of the IDM (not very stable)
  - 4: Macroscopic version of the OVM (not very stable, obsolete)
  - 5: Black Scholes equation (yes, MASTER can also do some finance math for fixed-strike-time derivatives!). Then, the model parameters v0->volatility, Tr=risk-free interest rate, tau0=stock price, xmax=stock price interval, choice_IC for the payment profile at strike (the simulation goes backwards) and choice_BC for kind of derivative (see demo project)q
  - 6: GKT model with "resignation effect"
  - ...
  - 8: OVM/FVDM/GKT with indifference zone a la "3-phase theory"
  - 50: Fokker-Planck equation


Basically, only the model choices 0 and 2 are useful. Notice that, with the GKT model, `choice_model`=0, one can also "mimick" first-order models by strongly reducing the relaxation time tau thus effectively fixing the flow/speed to the steady-state dictated by density, and also reducing the anticipation (but not to zero because, then, macroscopic second-order models are always flow (not numerically) unstable (see the demo projects)

* Choice of the explicit numerical integration method. `choice_method`==...

  - 0: Second-order McCormack: Needed for local models
  - 1: First-order upwind: For nonlocal models, the "downwind" information is ensured by the anticipative nonlocalities, so the upwind method is enough (and often better/more stable)

* Besides two test flags, there is also an optional line for `dx_rmp` which is only used if `choice_rmp==1` and therefore obsolete.



## Optional or conditional-mandatory project files
------------------------------------------------

### `proj.v0` or `proj.Tr`

Sets a space dependency of the parameters desired speed v0 or desired time gap in following Tr thereby overriding posible settings in the `proj.inp` file. This serves to implement flow-conservative bottlenecks. For example, a .v0 file can look like this:

```
%  x(m)    v0 (m/s)
%------------------
   4600    33
   5000    28  
```

Lines beginning with `%` are comment lines.
The data lines say that,
  * for x<4600 m, the desired speed is 33 m/s
  * for x>5000 m, it is 28 m/s
  * in between, it changes linearly from 33 m/s to 28 m/s
Particularly, this means that jumps of the values are implemented by two lines with the same location and the old and new v0 value

### `proj.A`

has only two entries pos_rel and width_rel, which, together with the input parameters A0 and dA define the squared GKT speed variation coefficient curve
A(rho)=A0 + dA * ( tanh(arg)+1), arg=(rho/rhomax - pos_rel)/width_rel


### `proj.IC`

has the columns x[m], rho[veh/km], and Q[veh/h]. Needed for `choice_init`=3 or =4. Interpolation/extrapolation as for the .v0 and .Tr files. 

### `proj.BCl` and `proj.BCr`

They have the columns time [s], rho [veh/km] and flow [veh/h]. Interpolation/extrapolation as for the .v0 and .Tr files.

Watch out for some attempted detection in the program whether SI units have been used at input. Better use veh/km and veh/h throughout for the .IC, .BCl, and .BCr files

### `proj.rmps`

Has the columns center merge/diverge position [m], and length Lrmp [m]. For each data line n, a file `proj.rmp<n>` specifying the flows is needed

### `proj.rmp<n>`, n=1,2,...

Has the two columns time [s] and ramp flow [veh/h]. Inter- and extrapolation as in the .v0 file.

### `proj.KKL`

model parameters for the Kerner-Konhaeuser (-Lee) model. (`choice_model`=2, self-explaining

### `proj.IDM`

Additional parameters for macroscopoic IDM (`choice_model`=3). Only the parameters b, s0, and delta are obtained since v0, a=v0/tau, and T=Tr are already defined by the `inp` file 

### `proj.3PHASE`

Self-explaining. For `choice_model=8`. Has its own `choice_variant` entry selecting which model should be brought into a 3-phase version. Overseeds the `.OVM` file

### `proj.FPE`

For the Fokker-Planck equation (`choice_model=50`)



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


