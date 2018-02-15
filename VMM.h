#ifndef VMM_H
#define VMM_H

#include "master.const"


/** second-order model from Velasco and Marques, submitted to
PRE (aug-05)
choice_model=21
*/

class VMM{

 public:

  VMM();
  VMM(double dx, double xWidth, double vehLen, int choice_method);

  void get_modelparams(const char projName[]);

  void calc_rhs(bool downwind_diff,
              const double rho[],const double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[]);

 private:

  void calc_eq();
  char fName[256];
  double veqtab[NRHO+1];
  double dQdx[NXMAX+1]; //dQ/dx
  double dx;
  int nx;

  double v0;
  double tau0;
  double alpha; //tau=tau0*(1+alpha*(rho/rhomax)^2)
  double theta0; //theta=theta0*(1-rho/rhomax)
  double rhocrel; //Ve=v0*(1+exp(rho/rhomax-rhocrel)/drhorel)^{-1}-3.72e-6)  
  double drhorel; 
  double rhomax;

};

#endif // VMM_H


