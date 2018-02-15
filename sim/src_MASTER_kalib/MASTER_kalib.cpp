

// c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

//alternatively there is <cstdio> which declares
//everything in namespace std
//but the explicit "using namespace std;" puts
//everything in global namespace


// c++ 
#include <iostream>
using namespace std;

// own
#include "general.h"

#include "Statistics.h"
#include "RandomUtils.h" // contains, e.g.,  myRand()
#include "InOut.h"
//#include "TemplateClass.h"

// constants

static const int NDATA_MAX=5000;// max. number of data points
static const int MAXSTR=500;// max. number of data points

// Example helper function w/o .h files

double objFunOneDet(int ndata, double rho[NDATA_MAX], double V[NDATA_MAX], 
		    double rhoref[NDATA_MAX], double Vref[NDATA_MAX]){
   double mu=1.;
   if(ndata>NDATA_MAX){
    cerr<<"Error: number of data points ndata>NDA^TA_MAX"<<endl;
    exit(-1);
   }
   double sum=0;
   for (int i=0; i<ndata; i++){
     sum += mu*pow(1-V[i]/Vref[i], 2) + (1-mu)*pow(1-rho[i]/rhoref[i], 2);
   }
   return sum/ndata;
}


 
//#####################################################
//#####################################################

 //#####################################################
 // main
 //#####################################################

int main(int argc, char* argv[]) {

  static const int NDATA_MAX=5000;// max. number of data points
  static const int N2D_MAX=100;
  static const int NX_MAX=10;
  static const int MAXSTR=200; // max string length

  double incr_T=0.1;
  double Tmin=1.0;
  double Tmax=2.0;
  int nT=(int)((Tmax-Tmin)/incr_T+1.5);

  double incr_tau=2;
  double taumin=10;
  double taumax=30;
  int ntau=(int)((taumax-taumin)/incr_tau+1.5);

  int incr_x=2000;
  int xmin=10000;
  int xmax=16000;
  int nx=(xmax-xmin)/incr_x+1;
  int xvals[N2D_MAX];
  for (int ix=0; ix<nx; ix++){xvals[ix]=xmin+incr_x*ix;}

  double rhoref[NX_MAX][NDATA_MAX];
  double Vref[NX_MAX][NDATA_MAX];
  double rhodata[NDATA_MAX];
  double Vdata[NDATA_MAX];
  int ndata;
  char projRefName[MAXSTR];
  char projName[MAXSTR];
  char detFileName[MAXSTR];

//#####################################################
  double Tref=1.5;
  double tauref=16;
//#####################################################

  sprintf(projRefName,"kalib_T%.1f_tau%.0f",Tref,tauref);

  InOut inout;

  for (int ix=0; ix<nx; ix++){
    cout <<"xvals[ix]="<<xvals[ix]<<endl;
    sprintf(detFileName,"%s.x%i", projRefName, xvals[ix]);
    inout.get_col(detFileName,3,ndata, rhoref[ix]);
    inout.get_col(detFileName,4,ndata,Vref[ix]);
  }

  cout <<"HIER"<<endl;

 //#####################################################
 // Calculation of objective-fun landscape
 //#####################################################

  double objFun[N2D_MAX][N2D_MAX];

  for (int iT=0; iT<nT; iT++){ 
    for (int itau=0; itau<ntau; itau++){ 
      double T=Tmin+iT*incr_T;
      double tau=taumin+itau*incr_tau;
      sprintf(projName,"kalib_T%.1f_tau%.0f",T,tau);
      cout <<"projName="<<projName<<endl;

      objFun[iT][itau]=0;

      for (int ix=0; ix<nx; ix++){
	sprintf(detFileName,"%s.x%i", projName, xvals[ix]);
	inout.get_col(detFileName,3,ndata,rhodata);
	inout.get_col(detFileName,4,ndata,Vdata);
	objFun[iT][itau]+= objFunOneDet(ndata, rhodata, Vdata, rhoref[ix], Vref[ix]);

      }
    }
  }

 //#####################################################
 // output
 //#####################################################


  char outfileName[MAXSTR];
  FILE  *outfile;                  
  sprintf(outfileName,"kalib_objFun.dat");

  cout <<"writing objective function landscape in "<<outfileName<<" ..."<<endl;
  outfile = fopen(outfileName,"w");
  fprintf(outfile, "#T\ttau\tobjFun\n");

  for (int iT=0; iT<nT; iT++){ 
    for (int itau=0; itau<ntau; itau++){ 
      double T=Tmin+iT*incr_T;
      double tau=taumin+itau*incr_tau;
      fprintf(outfile,  "%.2f\t %.2f\t %.6f\n", T,tau,objFun[iT][itau]);
    }
    fprintf(outfile, "\n");        // newline mean new "sweep" for gnuplot
  }
  fclose(outfile);


 cout<<"\nMASTER_kalib main call finished ... \n\n";
 return(0);
}

// #########################################################







