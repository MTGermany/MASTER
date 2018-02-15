#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// c++ 
#include <iostream>
#include <fstream>

using namespace std;

// own

#include "VMM.h"
#include "general.h"
#include "InOut.h"
//#include "../statistics/Statistics.h"
VMM::VMM(){}

VMM::VMM(double dx, double xWidth, double vehLen, int choice_method){
  if(choice_method==1){
    error(" VMM is parabolic; use McCormack instead of upwind");
  }

  this->dx=dx;
  nx   =(int)(xWidth/dx);
  rhomax = 1./vehLen; 

  if(false){
    cout <<"VMM Cstr: dx="<<dx<<" nx="<<nx<<endl;
    //exit(0);
  }
}


void VMM::get_modelparams(const char projName[]){

  sprintf(fName,"%s.VMM",projName);

  FILE *fp;
  InOut inout;
  fp=fopen(fName,"r");
  if(fp){
    filecheck(fp,fName);

    inout.getvar(fp,&v0); 
    inout.getvar(fp,&tau0); 
    inout.getvar(fp,&alpha); 
    inout.getvar(fp,&theta0); 
    inout.getvar(fp,&rhocrel); 
    inout.getvar(fp,&drhorel); 
    fclose(fp);
  }
  if(false){
    cout <<"FPE get_modelparams: v0="<<v0<<" drhorel="<<drhorel
         <<endl;
    //exit(0);
  }
  calc_eq();

}

void VMM::calc_eq(){
  for(int ir=0; ir<=NRHO; ir++){
    double rhorel=ir/NRHO;
    veqtab[ir]=v0*(1/(1+exp(rhorel-rhocrel)/drhorel)-3.72e-6);
  }
}



void VMM::calc_rhs(bool downwind_diff,
              const double rho[],const double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[]){

  int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                        // ishift=1 for downwind diffrences

  for (int i=0; i<=nx; i++){
    dQdx[0]=0;
    dQdx[nx]=0; // no diffusion at the boundaries 
    for (i=1; i<=nx-1; i++){
       dQdx[i] = ( Q[i+ishift]-Q[i+ishift-1])/dx;
    }
  }

  for (int i=0; i<=nx; i++){
      // v0fac = ve(rho)/ve(0) -> v0factab calculated in equil.cc
    double rholoc=(rho[i]>0.0001) ? rho[i]  : 0.0001;

    double Ve = intp(veqtab, NRHO+1, rholoc,0,rhomax); 
    double Qe=rholoc*Ve;
    double thetae=theta0*(1-rholoc/rhomax);
    double taue=tau0/(1+alpha*SQR(rholoc/rhomax));
    double eta=rholoc*taue*Qe;
    F1[i] = Q[i];
    F2[i] = (SQR(Q[i])- eta*dQdx[i])/rholoc + rholoc*thetae;
    S1[i] = 0;
    S2[i] = (rholoc*Ve - Q[i]) / taue; 
    if(false){cout <<"VMM.calc_rhs: nx="<<nx
		  <<" x="<<(i*dx)
		  <<" rho="<<rholoc
		  <<" Q="<<Q[i]
                  <<" F1="<<F1[i]
                  <<" S1="<<S1[i]
                  <<" F2="<<F2[i]
                  <<" S2="<<S2[i]
		  <<endl;
    }
  }

}



