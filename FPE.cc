
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// c++ 
#include <iostream>
#include <fstream>


using namespace std;

// own

#include "FPE.h"
#include "general.h"
#include "InOut.h"


FPE::FPE(){}

FPE::FPE(double dx, double xWidth,int choice_method){
  if(choice_method==1){
    error(" FPE is parabolic; use McCormack instead of upwind");
  }

  this->dx=dx;
  this->xWidth=xWidth;
  nx   =(int)(xWidth/dx);
  if(false){
    cout <<"FPE Cstr: A="<<A<<" B="<<B<<" nx="<<nx<<endl;
    //exit(0);
  }
}

void FPE::get_modelparams(const char projName[]){

  sprintf(fName,"%s.FPE",projName);

  FILE *fp;
  InOut inout;
  fp=fopen(fName,"r");
  if(fp){
    filecheck(fp,fName);

    inout.getvar(fp,&xmin); // xWidth from .inp file (called xmax there)! 

    inout.getvar(fp,&A);    // drift coefficient
    inout.getvar(fp,&B);    // second Kramers-Moyal coeffixient (2*Diff)
    fclose(fp);
  }
  else{cerr<<" file "<<fName<<" does not exist"<<endl; exit(-1);}
}



void FPE::calc_rhs(const double f[],
              double F1[],
              double S1[],
              bool downwind_diff){

    int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                        // ishift=1 for downwind diffrences

    for (int i=0; i<=nx; i++){
      F1[i]= A* f[i] - 0.5 * B * ( f[i+ishift]-f[i+ishift-1])/dx;
      S1[i]=0;
    }
}



