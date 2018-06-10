#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// c++ 
#include <iostream>
#include <fstream>

using namespace std;

// own

#include "SGM.h"
#include "general.h"
#include "InOut.h"
//#include "../statistics/Statistics.h"
SGM::SGM(){}

SGM::SGM(double dx, double xWidth, double vehLen, int choice_method){
  if(choice_method==1){
    error(" SGM is parabolic; use McCormack instead of upwind");
  }

  this->dx=dx;
  nx   =(int)(xWidth/dx);
  rhomax = 1./vehLen; 

  if(false){
    cout <<"SGM Cstr: dx="<<dx<<" nx="<<nx<<endl;
    //exit(0);
  }
}


void SGM::get_modelparams(const char projName[]){

  sprintf(fName,"%s.SGM",projName);

  FILE *fp;
  InOut inout;
  fp=fopen(fName,"r");
  if(fp){
    filecheck(fp,fName);

    inout.getvar(fp,&v0); 
    inout.getvar(fp,&rhomax); 
    inout.getvar(fp,&rhocrel); 
    inout.getvar(fp,&drhorel); 
    inout.getvar(fp,&c0); 
    inout.getvar(fp,&tau); 
    fclose(fp);
  }

  calc_eq(veqtab); // !! veqtab local and global parallel!
  if(false){
    cout <<"SGE get_modelparams: v0="<<v0<<" drhorel="<<drhorel
         <<endl;
    //exit(0);
  }
}

void SGM::calc_eq(double veqtab[]){
  for(int ir=0; ir<=NRHO; ir++){
    double rhorel=((double)ir)/NRHO;
    veqtab[ir]=v0*(1/(1+exp((rhorel-rhocrel)/drhorel))-3.72e-6);
    cout<<"SGM: ir="<<ir<<" rhorel="<<rhorel<<"  veqtab="<< veqtab[ir]<<endl;
  }
}



void SGM::calc_rhs(bool downwind_diff,
              const double rho[], const double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[]){

  int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                        // ishift=1 for downwind diffrences
  double v[NXMAX+1];

  for (int i=0; i<=nx; i++){
      v[i] = Q[i]/rho[i];
  }

     // addl ad-hoc diffus!

  double flux_diffus[NXMAX+1];  // diffusion comes into the flux F2
  double mu=166.7; //as KKL model 166.7
  flux_diffus[0]=flux_diffus[nx]=0; // no diffusion at the boundaries 
  for (int i=1; i<=nx-1; i++){ 
     flux_diffus[i] = - mu*( v[i+ishift]-v[i+ishift-1])/dx;
  }

  // calculations

  for (int i=0; i<=nx; i++){

    double Ve = intp(veqtab, NRHO+1, rho[i],0,rhomax); 
    double Qe=rho[i]*Ve;
    double P=-c0*Q[i+ishift];


    // output

    F1[i] = Q[i];
    F2[i] = SQR(Q[i])/rho[i] - P + flux_diffus[i];
    S1[i] = 0;
    S2[i] = (Qe - Q[i]) / tau; 
    if(false){
    //if(i<5){
      cout <<"SGM.calc_rhs: "
		  <<" x="<<(i*dx)
		  <<" rho="<<rho[i]
		  <<" Q="<<Q[i]
		 <<" Qe="<<Qe
                  <<" F1="<<F1[i]
                  <<" S1="<<S1[i]
                  <<" F2="<<F2[i]
                  <<" S2="<<S2[i]
		  <<endl;
    }
  }

}



