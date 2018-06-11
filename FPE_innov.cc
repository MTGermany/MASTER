#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// c++ 
#include <iostream>
#include <fstream>


using namespace std;


// own

#include "FPE_innov.h"
#include "general.h"
#include "InOut.h"
// #include "Statistics.h"
FPE_innov::FPE_innov(){}

FPE_innov::FPE_innov(double dx, double xWidth,int choice_method){
  if(choice_method==1){
    error(" FPE_innov is parabolic; use McCormack instead of upwind");
  }

  this->dx=dx;
  this->xWidth=xWidth;
  nx   =(int)(xWidth/dx);
  tab_i=0;

  for (int i=0; i<=nx; i++){  // lookup table for x values
    x[i] = xmin + xWidth*i/nx; 
  }

  if(false){
    cout <<"FPE_innov Cstr: dx="<<dx<<" nx="<<nx<<endl;
    //exit(0);
  }
}

void FPE_innov::write_statProperties(const char projName[]){
  sprintf(fName,"%s.stat",projName);
  FILE *fp;
  InOut inout;
  fp=fopen(fName,"w");
  if(fp){
    inout.write_array (fName, tab_i-1,
		       tab_time,tab_erw_x,tab_var_x,tab_gamma_x,
		       "Time","<x>","Var(x)","Skewness(x)");
  }
  else cerr <<"Error: Could not write to "<<fName<<endl;
}

void FPE_innov::save_statProperties(const double f[], double time){

  // calculate skewness

  double gamma_x=0;
  for (int i=0; i<=nx; i++){
    gamma_x += SQR(x[i])*x[i]*f[i];
  }
  gamma_x /= norm;
  gamma_x -= (3*var_x+SQR(erw_x))*erw_x;

  tab_time[tab_i]=time;
  tab_erw_x[tab_i]=erw_x;
  tab_var_x[tab_i]=var_x;
  tab_gamma_x[tab_i]=gamma_x/pow(var_x, 1.5);

  tab_i++;
}


void FPE_innov::get_modelparams(const char projName[]){

  sprintf(fName,"%s.FPE",projName);

  FILE *fp;
  InOut inout;
  fp=fopen(fName,"r");
  if(fp){
    filecheck(fp,fName);

    inout.getvar(fp,&xmin); // xWidth from .inp file (called xmax there)! 

    inout.getvar(fp,&A);    // drift coefficient
    //A=0;
    inout.getvar(fp,&D);    // Diff, 0.5*second Kramers-Moyal coeffixient B

    inout.getvar(fp,&lambda);  //innovation rate
    inout.getvar(fp,&c0);  //if !=0, then multiplicative Diffusion D=c0*Var(x)
    inout.getvar(fp,&c1);  //if !=0, then multiplicative Diffusion D=c1*erw(x)
    fclose(fp);
  }
  else{cerr<<" file "<<fName<<" does not exist"<<endl; exit(-1);}

  if(false){
    cout <<"FPE get_modelparams: A="<<A<<" D="<<D<<" nx="<<nx
	 <<" lambda="<<lambda
	 <<" c0="<<c0
	 <<" c1="<<c1
         <<endl;
    //exit(0);
  }

}


void FPE_innov::normalize(double f[]){

  norm=0; 

  for (int i=0; i<=nx; i++){
    norm += f[i];
  }
  for (int i=0; i<=nx; i++){
    f[i] /= norm*dx;
  }
}


void FPE_innov::calc_rhs(const double f[],
              double F1[],
              double S1[],
              bool downwind_diff){

  // (1) calculate erw(x) and var(x)  (stat methods not applicable!)
  //      for effectivity: only necessary if (!)downwind_diff =>later

  if(true){ 
    erw_x=0;
    var_x=0;
    norm=0; // f should be normalized, but you never knows ...; also test

    for (int i=0; i<=nx; i++){

      norm += f[i];
      erw_x += x[i]*f[i];
      var_x += SQR(x[i])*f[i];
    }
    erw_x /= norm;
    var_x /= norm;
    var_x -= SQR(erw_x);

    if(false){
      cout <<" FPE_innov.calc_rhs: erw_x="<<erw_x<<" var_x="<<var_x
	   <<endl;
    }

  }

  // (2) Possibly re-normalize f exactly to 1 (then no attribute "const f")
  //     => check later

  // (3) calculate possibly multiplicative diffusion constant
  // if D constant, nothing need to be done

  if(fabs(c0)>1.e-6){
    D=c0*var_x;
    //cout <<" c0!=0 => D=c0*var_x="<<D<<endl;
  }
  if(fabs(c1)>1.e-6){
    D=c1*erw_x;
    //cout <<" c1!=0 => D=c1*erw_x="<<D<<endl;
  }



  // (4) do actual update step

  int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                        // ishift=1 for downwind diffrences

  for (int i=0; i<=nx; i++){

    F1[i]= A* f[i] - D * ( f[i+ishift]-f[i+ishift-1])/dx; //A=0
    S1[i]=lambda*f[i]*(x[i]-erw_x);
  }

}



