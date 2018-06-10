

/*########################################################################
  MASTER
''Macroscopic Simulation of Traffic to Enable Road Predictions''

(c) Martin Treiber, Feb 2000
#########################################################################

Traffic model: 
-----------------
     Gaskinetic-based traffic model (GKT model):
     macroscopic  equations for the macroscopic
     density rho(x,t), and the macroscopic velocity v(x,t).
     Effective-one-lane treatment, but on-ramps and off-ramps included.

     For the GKT model equations, see 
     PRL 81, 3042 (1998) and PRE 59, 239 (1999)

Integration schemes:
--------------------
     Mc-Cormack  (second consistency order in t and x)   
     Upwind      (first consistency order in t, second in x)

   - fixed grid with fixed time steps. 
   - Conservation of vehicles is automatically satisfied by using the
     flux representations of the traffic equations in the form
       d_t rho(x,t) = -d_x Q  + S1,
       d_t   Q(x,t) = -d_x F2 + S2. .h
     (cf. R.J.LeVeque, "Numerical Methods for Conservation Laws", 
     Birkhaeuser)
   ->Dynamical fields are the density rho and traffic flow Q=rho*v
   - For the GKT model, the upwind integration scheme
     is best. Mc-Cormack or other second-order schemes
     best for models without nonlocalities
     (Kerner-Konhaeuser-Lee etc),  

Further Info:
-----
   see file READMEmaster

######################################################################### */

/* c */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// c++ 
#include <iostream>
#include <fstream>
using namespace std;


// own
#include "general.h"
#include "master.const"
#include "master.h"

   using namespace std;
//************************************************************************
//************************************************************************

int main(int argc, char *argv[])

{
double rho[NXMAX+1];  // Density rho[0]...rho[nx] at the Actual time step 
double   Q[NXMAX+1];  // Flow rho*v 
double acc[NXMAX+1];  // Effective local acceleration = S2/Q


// Following arrays contain only output information
// at fixed locations ("simulated loop data")

double rho_xt[NXCUTMAX+1][NTLOOPMAX+1];
double   Q_xt[NXCUTMAX+1][NTLOOPMAX+1]; 
double   a_xt[NXCUTMAX+1][NTLOOPMAX+1]; 


// ##############################################################
// Input
// ##############################################################

if (argc!=2){
  cout <<"\nCalling sequence: master <simulation_set>," << endl
       << "where <simulation_set> is entered w/o extensions.";
  exit(-1);
}
else sprintf(projName,"%s",argv[1]);


// Read in all necessary parameters

input(projName);

// ##############################################################
// Initialization
// ##############################################################

// initialize models apart from GKT

 cout <<"projName="<<projName<<endl;
 
if(choice_model==8){
  mac3phases=Mac3phases(projName,xmax,dx);
}
 
if(choice_model==50){
  fpe=FPE(dx,xmax,choice_method);
  fpe.get_modelparams(projName);
}
if(choice_model==51){
  fpe_innov=FPE_innov(dx,xmax,choice_method);
  fpe_innov.get_modelparams(projName);
  fpe_innov.normalize(rho);
}
if(choice_model==21){
  vmm=VMM(dx,xmax,1./rhomax, choice_method);
  vmm.get_modelparams(projName);
}
if(choice_model==10){
  sgm=SGM(dx,xmax,1./rhomax, choice_method);
  sgm.get_modelparams(projName);
}

 
// Calculate tables of equilibrium velocity and flow, and of
// the interaction function B(delta_V) of the GKT model

calc_tables(projName);

int    nt     = (int) (tmax/dt);            // Iteration in time: it=0..nt
int    nx     = (int) (xmax/dx);            // Iteration in space: it=0..nx
int    ndtout = ((dtout/dt)>=1) 
              ? (int) (dtout/dt) : 1;       // Every ndtout'th timestep saved

// set initial conditions

initialize(rho,Q);


// save initial values

write_results (projName,0,rho,Q,acc);     

if (choice_outp ==1) {
  if (fabs(time_slices[0]) < 0.5*dt) 
    write_cross_section (1, 0, 0, rho,Q,acc);
  for (int ixcuts=0; ixcuts<=n_xcuts; ixcuts++)
  {
    rho_xt[ixcuts][0] = rho[(int) ((x_xcuts[ixcuts]-xmin)/dx)];
    Q_xt  [ixcuts][0] = Q  [(int) ((x_xcuts[ixcuts]-xmin)/dx)];
    a_xt  [ixcuts][0] = 0;
  }
}




// ##############################################################
// Main loop over the simulation time
// ##############################################################

cout << endl;

for (it=1; it<=nt; it++)
{

  // The actual simulation

  timestep (it,rho,Q,acc);

  // 3D output to file and control output to stdout

  if(it % ndtout==0) 
  {
    write_results (projName,it,rho,Q,acc);
    cout << "t = " << it*dt;

    // test of conservation of number of vehicles if periodic BC

    if (choice_BC==0)
    {
      double n_veh=0;    
      for (int ix=1; ix<=nx-1; ix++) n_veh+=rho[ix]*dx; // Length (n-1)*dx;
      cout << " n_veh = " << n_veh;
    }

    cout << endl;
  }

  // Optional additional output of time series (fixed x) and snapshots 
  // (fixed t). Assumes that simulation time step dt is always below 1 s

  if(choice_outp ==1)
  { 
    // write snapshot in a file

    for(int i=0; i<=n_slices; i++){
      if (fabs(time_slices[i] - it*dt) < 0.5*dt){
        write_cross_section (1, it*dt, 0, rho,Q,acc);
      }
    }

    // save data of time series for future writing in a file

    int itloop = (int)(it*dt/dtloop); 
    if (fabs( it*dt - itloop*dtloop) <= dt){ // it*dt always > itloop*dtloop 
      for (int ixcuts=0; ixcuts<=n_xcuts; ixcuts++){
        rho_xt[ixcuts][itloop] = rho[(int) ((x_xcuts[ixcuts]-xmin)/dx)];
        Q_xt  [ixcuts][itloop] = Q  [(int) ((x_xcuts[ixcuts]-xmin)/dx)];
        a_xt  [ixcuts][itloop] = acc[(int) ((x_xcuts[ixcuts]-xmin)/dx)];
      } 
      if(choice_model==51){
	fpe_innov.save_statProperties(rho,it*dt);
      }
    }

  }

}

// ##############################################################
// End simulation
// ##############################################################


// write optional time series data to file

if(choice_outp ==1){

  double rho_t[NTLOOPMAX+1]; 
  double   Q_t[NTLOOPMAX+1]; 
  double   a_t[NTLOOPMAX+1]; 

  for (int ix=0; ix<=n_xcuts; ix++)
  {
    for (int it=0; it<=(int)(nt*dt/dtloop); it++) 
    {
      rho_t[it] =  rho_xt[ix][it];
      Q_t  [it] =  Q_xt  [ix][it];
      a_t  [it] =  a_xt  [ix][it];
    }
    write_cross_section (0, 0, (int) x_xcuts[ix], rho_t,Q_t,a_t);
  }
}

// write optional statistical properties if FPE_innov model
if(choice_model==51){
  fpe_innov.write_statProperties(projName);
}


return(0);
}

#include "general.cc"
#include "io.cc"
#include "equil.cc"
#include "initialize.cc"
#include "bc.cc"
#include "timestep.cc"
#include "pde.cc"
#include "calc_rhs.cc"
#include "calc_rhsGKT.cc"
#include "KKL.cc"
#include "macOVM.cc"
#include "macIDM.cc"
#include "trucks.cc"







