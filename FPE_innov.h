#ifndef FPE_INNOV_H
#define FPE_INNOV_H

#include "master.const"


class FPE_innov{

 public:

  FPE_innov();
  FPE_innov(double dx, double xWidth,int choice_method);

  void get_modelparams(const char projName[]);
  void write_statProperties(const char projName[]);
  void save_statProperties(const double f[], double time);
  void normalize(double f[]);

  void calc_rhs(const double f[],
              double F1[],
              double S1[],
		bool downwind_diff);

 private:


  char fName[256];

  double A;
  double D; // =B/2

  double lambda; //imitation rate
  double c0; //if c0 !=0, then multiplicative Diffusion D=c0*Var(x)
  double c1; //if c1 !=0, then multiplicative Diffusion D=c1*erw(x)

  double xmin;
  double xWidth;
  double dx;
  int nx;

  double erw_x;  // expectation value int dx x f(x)
  double var_x;  // variance

  int tab_i;
  double tab_time[NTMAX];  // tables for further saving
  double tab_erw_x[NTMAX];
  double tab_var_x[NTMAX];
  double tab_gamma_x[NTMAX];  // skewness <(x-erw_x)^3>/var_x^{3/2}
  double norm; // f should be normalized, but you never knows ...; also test
  double x[NXMAX]; // lookup table 
};


#endif // FPE_INNOV_H

