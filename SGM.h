#ifndef SGM_H
#define SGM_H

#include "master.const"


/** second-order model from Velasco and Marques, submitted to
PRE (aug-05)
choice_model=21
*/

class SGM{

 public:

  SGM();
  SGM(double dx, double xWidth, double vehLen, int choice_method);

  void get_modelparams(const char projName[]);

  void calc_rhs(bool downwind_diff,
              const double rho[], const double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[]);
  void calc_eq(double veqtab[]);

 private:

  char fName[256];
  double veqtab[NRHO+1];
  double dQdx[NXMAX+1]; //dQ/dx
  double dx;
  int nx;

  double v0;
  double rhomax;
  double rhocrel; //Ve=v0*(1+exp(rho/rhomax-rhocrel)/drhorel)^{-1}-3.72e-6)  
  double drhorel; 
  double c0;  // additional V d_x V term: (V-c0) d_x V
  double tau; 

};

#endif // SGM_H


