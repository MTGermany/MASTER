#ifndef Mac3phases_H
#define Mac3phases_H


class Mac3phases{

 public:

  Mac3phases();
  Mac3phases(const char projName[],
	     double xWidth, double dx); 
  
  void get_modelparams(const char fname[]);

  void calc_rhs (
              double rho[],double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[],  int choice_BC,
		 int show_calc_rhs);

  double get_v0(){return v0;}

  double rhoFree_Q(double Q);
  double rhoCong_Q(double Q);
  double rho_v(double v);
  double ve_rho(double rho);

 
 private:
  int counter;
  bool intp_err;
  void calc_tables();
  void write_tables();
  double intp(const double tab[], int n,double x, 
	      double xmin, double xmax);

  void shift_in_x (double antic_factor,  int choice_BC,
		 const double v[], const double f[], double fshifted[]);
  
  double vopt(double seff,double v);
  
  //void add_ramps(const double rho[], const double Q[], double S1[], double S2[]);
  //void add_rmpGKT ( const double rho[], const double Q[],
  //		    double S1[], double S2[],
  //		    double x_rmp, double dx_rmp, double flow_rmp);

  char projName[256];
  
  double Atab[NRHO+1];   // A = Dimensionless measure for vel variance
  double veqtabmin[NRHO+1]; // eq. velocity for Tmax
  double veqtabmax[NRHO+1]; // eq. velocity for Tmin
  double rho_vtab[NRHO+1]; // rho_e(v) ( inverse of veqtabmax)
  double rho_freetab[NRHO+1]; // free branch of rho_e(Q)
  double rho_congtab[NRHO+1]; // congested branch of rho_e(Q)

  double Arhomax;
  double rhoQmax;
  double Qmax;
  double Qtabmax;
  
  double xmax;
  double dx;
  int nx;

  int choice_variant;
  double rhomax;
  double v0;
  double tau;
  double Tmin;
  double Tmax;
  double s0;       // addtl. introduced distance (m) for veq=0
  double beta;    // addtl. introduced distance (m) for veq=0

    
  double A0,dA;
  double posA_rel,widthA_rel;

  double lveh;
  double seff;


};

#endif // Mac3phases_H
