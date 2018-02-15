
// 
// original KK equilibrium relation (Hermann-Kerner preprint. 1997):
// ve = v0/(1.+exp((rho/rhomax-rhoi/rhomax)/b)) - d;
// where  d= 1./(1.+exp((1. -rhoi/rhomax )/b))

// Ve(rho) relation of Lee/Lee/Kim, PRL81,1130 (1998):
// ve = v0 (1 - rho/rhomax) / (1 + E ( rho/rhomax) ^theta)

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

using namespace std;

// c++ 
#include <iostream>
#include <fstream>


using namespace std;

// own

#include "KKL.h"
#include "general.h"
#include "master.const"
#include "master.h"

KKLparams::KKLparams()
{
  c0 = 15.;           // kerner: 41.94/3.6 m/s
  mu = 166.7;        // Kerner: 204/3.6  m/s
  tau= 30.;          // Kerner: 5 s
  lee_E  = 100.; 
  lee_theta = 4;
  kerner_rhoi = 38.5;
  kerner_b = 0.042;
}


void KKLparams::getKKLparams(char namepar[])
{
  FILE * fp;                    
  char   param_file[MAXSTR];   
  sprintf(param_file,"%s.KKL",namepar);
  cout << "\nKKL model chosen:" << endl
       << "Reading params of KKL model from file "<< param_file
       << endl;
  fp = fopen(param_file,"r");
  filecheck(fp, param_file);

  getvar(fp,&c0);
  getvar(fp,&mu);
  getvar(fp,&tau);
  getvar(fp,&lee_E); 
  getvar(fp,&lee_theta);
  getvar(fp,&kerner_rhoi);
  getvar(fp,&kerner_b);

  fclose(fp);
}
