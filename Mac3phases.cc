
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// c++ 
#include <iostream>
#include <fstream>


using namespace std;

// own

#include "master.const"
#include "Mac3phases.h"
#include "general.h"
#include "InOut.h"


Mac3phases::Mac3phases(){}

Mac3phases::Mac3phases(const char projName1[],
		       double xWidth, double dx){

  counter=0;
  intp_err=false;
  cout <<"in Mac3phases file Cstr: projName="<<projName1<<endl;


  sprintf(projName,"%s",projName1);
  
  this->xmax=xWidth;
  this->dx=dx;
  nx   =(int)(xmax/dx);
  
  char fName[256];
  sprintf(fName,"%s.3PHASE",projName);

  get_modelparams(fName);
  calc_tables();
  write_tables();
  if(true){
    cout <<"Mac3phases file Cstr: nx="<<nx<<" widthA_rel="<<widthA_rel<<endl;
    //exit(0);
  }
}

void Mac3phases::get_modelparams(const char fName[]){


  FILE *fp;
  InOut inout;
  fp=fopen(fName,"r");
  if(fp){
    filecheck(fp,fName);

    inout.getvar(fp,&choice_variant); // {original OVfun, triang. OV, 3phase}

    inout.getvar(fp,&rhomax); // vehicle length 
    inout.getvar(fp,&s0);   // addtl. introduced distance (m) for veq=0
    inout.getvar(fp,&v0); // martin mar08: Now always des. vel
    inout.getvar(fp,&tau);
    inout.getvar(fp,&Tmin);
    inout.getvar(fp,&Tmax);
    inout.getvar(fp,&beta); // prefactor of v*dv/s(rho) term

    inout.getvar(fp,&A0);  
    inout.getvar(fp,&dA);  
    inout.getvar(fp,&posA_rel);  
    inout.getvar(fp,&widthA_rel);
 
    fclose(fp);
  }
  else{cerr<<" file "<<fName<<" does not exist"<<endl; exit(-1);}

  Arhomax = A0 + 2.*dA; 
  lveh=1./rhomax;

}

void Mac3phases::calc_tables(){
  cout <<" in Mac3phases.calc_tables() ..."<<endl;
   for(int ir=0; ir<NRHO+1; ir++){
      double rho = rhomax*ir/NRHO;
      double s=1./rho - 1./rhomax;
      double arg      = (rho/rhomax - posA_rel)/widthA_rel; 
      Atab  [ir]  = A0 + dA * ( tanh(arg)+1.);
      veqtabmax[ir]=max(0., min(v0, (s-s0)/Tmin));
      veqtabmin[ir]=max(0., min(v0, (s-s0)/Tmax));
   }

   // maximum of flow
   
   rhoQmax=1/(lveh+s0+v0*Tmin);
   Qmax=v0*rhoQmax;
  
   // tables of rho(v), rhoFree (Q), rhoCong(Q)
   for (int iv=0; iv<NRHO+1; iv++){
     rho_vtab [iv] = 1/(lveh+s0+v0*iv/NRHO*Tmin);
   }
   double Qtabmax   = QFACTAB*Qmax;     // max. Q value of tables (>Qmax!)

   for (int iQ=0; iQ<NRHO+1; iQ++){
      double Qloc=Qtabmax*iQ/NRHO;
     rho_freetab[iQ] =(Qloc<=Qmax) ? Qloc/v0 : rhoQmax;
     rho_congtab[iQ] =(Qloc<=Qmax) ? (1-Qloc*Tmin)/(lveh+s0):rhoQmax;
   }
}


double Mac3phases::ve_rho(double rho){
  return intp(veqtabmax,NRHO+1,rho,0,rhomax);
}

double Mac3phases::rho_v(double v){
  //return intp(rho_vtab,NRHO+1,v,0,v0);
  return 1/(lveh+s0+v*Tmin);
}

double Mac3phases::rhoFree_Q(double Q){
  return intp(rho_freetab,NRHO+1,Q,0,Qtabmax);
}

double Mac3phases::rhoCong_Q(double Q){
  return intp(rho_congtab,NRHO+1,Q,0,Qtabmax);
}


//  function  intp GKT  
//#######################################################

double Mac3phases::intp(const double tab[], int n,double x, 
	     double xmin, double xmax)
   /* intp interpolates the array tab with n 
      equidistant points from 0 to n-1
      in [xmin,xmax] at the location x; an error message is produced on
      attempting extrapolation */
{
  double intp_value;
  double ir   = (n-1)*(x-xmin)/(xmax-xmin);
  int    i    = (int) ir;
  double rest = ir-i;
  if ((i>=0) && (i<n-1))  intp_value =  (1-rest) * tab[i] + rest*tab[i+1];
  else if (i==n-1) intp_value = tab[n-1];
  else {
    cout << "Mac3phases.intp: x="<<x<<" xmin="<<xmin<<" xmax="<<xmax
	 <<" n="<<n<<" index i = "<<i<<" (ir="<<ir<<") out of range\n";
     intp_value = tab[n-1];
     exit(-1);
     intp_err=true;
  }
  return(intp_value);
}


//#######################################################

void Mac3phases::write_tables(){
  char tabName[256];
  sprintf(tabName,"%s.tab1",projName);
  FILE *fptab;
  fptab=fopen(tabName,"w");
  cout <<"tabName="<<tabName<<endl;
  
  fprintf(fptab, "# rho(1/km)\t    A\t veqmax(km/h)\t Qeqmax(1/h)\t veqmin(km/h)\n");

  for(int ir=0; ir<NRHO+1; ir++) {
    double rho              = rhomax*ir/NRHO;
    fprintf(fptab, "%.1f\t  %.4f\t %.2f\t \t %.2f \t %.2f  \n",
       1000.*rho, Atab[ir], 
	    3.6*veqtabmax[ir], 3600.*rho*veqtabmax[ir], 3.6*veqtabmin[ir] );
    cout <<"rho="<<rho<<" rhomax="<<rhomax<<endl;
  }
  fclose(fptab);

}



void Mac3phases::calc_rhs (
			   double rho[],double Q[], 
			   double F1[], double F2[], 
			   double S1[], double S2[], int choice_BC,
			   int debug){

  counter++;

  double v[NXMAX+1];                         // velocity
  double rhodelta[NXMAX+1];                  // nonlocal density
  double vdelta[NXMAX+1];                    // nonlocal velocity


  /* Calculates the right-hand sides (fluxes Fi and sources Si)
     of macroscopic traffic models of the form
              ________________________________________________
              |  d_t rho(x,t)  = -d_x F1 + S1 (ramp flows)    |
              |  d_t Q(x,t)    = -d_x F2 + S2                 |
              ------------------------------------------------

     Input: 
       - choice_model: defines the traffic model; choice_model==0: GKT model
       - downwind_diff: whether downwind differences should be used
         for calculating derivatives (only relevant if fluxes with gradients)
       - The fields rho,Q. They will be modified only, if 
         something unphysical happens
       - switch "debug" to show some debugging information

     Output: 
       - The right-hand sides of the PDE:
         fluxes F1,F2, and sources S1,S2 of the model.
       - The sources include ramps: S1 = ramp flow/merging length 
                                    S2+= S1*(velocity when merging) 
       - Flow-conserving bottlenecks realized by gradients of the
         parameters T (=Tr) or v0, given by the (global) arrays
         Tr_loc and v0_loc
   */


  if( debug){
    printf("in calc_rhs; ");
  }


  // Check input (rho>0, Q>=0) and calculate V

  //bool stopped=false;
  
  for (int i=0; i<=nx; i++){
     if (rho[i]<TINY_VALUE)        rho[i] = TINY_VALUE;
     if (Q[i]<rho[i]*TINY_VALUE)   Q[i]   =rho[i]*TINY_VALUE;
     if (rho[i]>1/(lveh)){// !!: Keine Aenderung von rho, nur evtl. Q!
 
       //cout <<"Mac3phases.calc_rhs: stopped! x="<<i*dx<<" it="<<counter<<endl;

       //exit(-1);
     }
     v[i] = Q[i]/rho[i];
  }
  
  if  ( choice_variant==0){ // choice_variant=0: with A and equil nonlocalities

    // calculate advanced rho and v fields

    double antic_factor=1.0;
    shift_in_x(antic_factor, choice_BC, v, rho, rhodelta);
    shift_in_x (antic_factor, choice_BC, v, v,   vdelta);

    // calculate rhs of flow equation

    for (int i=0; i<=nx; i++){
       //double A     = intp(Atab, NRHO+1, rho[i], 0, rhomax);
       double A     = 0;
       double seff=max(2./(rho[i]+rhodelta[i])-lveh-s0, TINY_VALUE);

       F1[i] = Q[i];
       F2[i] = rho[i] * v[i]*v[i] * (1.+ A);

       double accOVM=(vopt(seff,v[i]) -v[i]) / tau;
       double dv=vdelta[i]-v[i];
       double accvdiff=beta*v[i]*dv/seff; 

       if(false){
          double b=2.;
          double exponent=2.;
          if(accvdiff<0){accvdiff=-pow(-accvdiff,exponent)/pow(b,exponent-1);}
       }
       double acc=max(-10*v0/tau, min(accOVM+accvdiff, v0/tau));
 
       S1[i]=0;
       S2[i] = rho[i]*acc;
      
    }

  } // choice_variant=0

  else if(choice_variant==1){
  }

  // add_ramps(rho,Q, S1, S2); // centrally in calc_rhs.cc
  
} // calc_rhs



double Mac3phases::vopt(double seff,double v){
  double small_val=1.e-6;
  double Tdyn=seff/max(v, small_val);
  return (Tdyn>Tmax)
      ? min(seff/Tmax, v0) : (Tdyn>Tmin)
      ? min(v+0.,v0) : (Tdyn>0)
      ? min(seff/Tmin, v0) : 0;
}



void Mac3phases:: shift_in_x (double antic_factor, int choice_BC,
		    const double v[], const double f[], double fshifted[]){

  //   Calculates, from the input array "f" representing f(x),
  //   an array "fshifted" representing f(x+dx_shift) 
  //   The shifted distance is s=antic_factor * (1/rhomax + Tr*v). 
  //   The extrapolation to the right depends on the BC:
  // !! changes here imply changes in function "boundary_vals" and vice versa

  double T=0.5*(Tmin+Tmax);

  if(antic_factor > 0){
    for(int i=0; i<nx+1;i++){
      double dx_shift  = antic_factor * ( 1/rhomax + s0+T*v[i] );
      int idx_shift = (int)(dx_shift/dx);
      double rest   = dx_shift/dx - idx_shift;

      if((i+idx_shift+1)<=nx){
        fshifted[i]    = (1-rest)*f[i+idx_shift] + rest*f[i+idx_shift+1];
      }
      // extrapolation depending on the BC: 0=periodic,1=free,2=Neumann, 3=Dirichlet

      else{ 
        if      (choice_BC==0) fshifted[i] 
            = (1-rest)*f[i+idx_shift-nx+1] + rest*f[i+idx_shift-nx+2];

        else if ((choice_BC==1) || (choice_BC==4)) fshifted[i] 
            = f[nx] + (dx_shift/dx - nx+i) * (f[nx]-f[nx-1]);

        else fshifted[i] 
	    = f[nx];  // const extrapolation in case of Neumann  or fixed BC
      }
    }
  }
  else{ // no shift
    for(int i=0; i<=nx;i++) fshifted[i]=f[i];
  }
}



// ramps centrally in calc_rhs.cc!


 
