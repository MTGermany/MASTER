

#ifndef MASTER_H
#define MASTER_H

#include "KKL.h"
#include "FPE.h"
#include "FPE_innov.h"
#include "VMM.h"
#include "SGM.h"
#include "Mac3phases.h"

//********************************************************************
//* global vars related to outer i/o and time loop
//********************************************************************

int    it;              // index of outer time loop
char   projName[MAXSTR]; // Name of the parameter set; only interactive action!

/********************************************************************/
/* Input parameters; Values in File "navst<npar>.inp" */
/********************************************************************/

/* model parameters */

double v0;             // Unique and constant desired velocity (m/s) 
double Tr;             // Reaction time (s) 
double tau0;           // Acceleration time at rho=0 (s) 
double A0;             // Acceleration fluctuation intens. at rho=0 
double rhomax;         // Max. density (the inverse of the average
                              // length of the vehicles (1/m) 
double antic_factor;   // distance of the "interaction point" 
double dA;             // amplitude of the variation of A

/* initial conditions */

double rho_init;       // Constant parts of the initial values 
double ampl_init;      // Ampl. of init. perturb. of rho (1/m) or Q (1/s) 
double w_init;         // width of (the envelope of) the perturbation (m)
double wavel_init;     // Wavelength of the sinusoidal perturbation (m) 
int    choice_init;    // {Q periodic*Gauss, shock, rho dipole-like}
double center_init;    // rel. center of init. perturbation (xmax)
int    choice_BC;      // {periodic,linear,Neumann,Dirichlet,Dutch}  
int    choice_rmp;     // {no ramp, one ramp, several ramps}
double dx_rmp;	       // merging length	

/* dynamic density sensitive velocity regulation */

int choice_dyn_v0;  //true if there is a file <>.v0_dyn
double x_v0[NBC];  //locations writen in <>.v0
int n_jmps_v0;  
double v0_x[NBC];  //velocities  "     "   " 
double delta_v0_x[NBC]; //  "      v-intervals by which v0 is lowered   
double thrsh_x[NBC];  //    "      rho - thresholds for v0-regulation

/* Range of integration and numerical parameters (xmin,tmin below) */

double tmax;           // Total time to be simulated (s) 
double dt;             // Largest possible discretisation in time (s) 
double dtmin;          // Smallest allowed discretisation in time 
double xmax;           // Length of the section to be simulated (m) 
double dx;             // Discretisation in space (m)  
double D1;             // Additional viscosity of the dens. eq. (m^2/s) 
double D2;             // Additional viscosity of the veloc. eq. (m^2/s) 


/* control of the output (what; in which form) */

double dtout;          // Discretisations for saving the results 
double dxout;          // (should be multiples of dt, dx) 
double dtloop;          // time interval of time series on output
int   choice_outp;     // {only equil+3D, +time and space cuts, +lin}  

/* Choice of numerical methods and model variants */

int choice_model;         // {v0,v0 with high-density-corr, vw_KK}
int choice_method;     // {McCormack pred_left, McC pred_right, upwind}
int choice_A;          // {A0, A0 + dA * inverse gauss function}  


/* Flags for tests  */

int test_timestep;
int test_shift;

// fields for the x-dependent (but t independent) model parameters

double v0_loc[NXMAX+1];
double Tr_loc[NXMAX+1];
// double antic_factor_loc[NXCUTMAX+1];


// fields for x and t dependent Tr (if flag_Txt is TRUE)

bool   flag_Txt;
int    n_Ts;      // how many x sections with T=Tr different from .inp
double x_Ts[NRMPMAX];             // center of these x sections
double dx_Ts[NRMPMAX];            // length of these x sections
double width_Ts[NRMPMAX];         // width of slopes of Tr
int    n_jumps_Ts[NRMPMAX];       // number of jumps in the sections
double times_Ts[NRMPMAX][NBC];    // times of jumps in the sections
double dT_Ts[NRMPMAX][NBC];       // differenc of Tr compared to .inp

// The coordinates for the start of the simulation
// ONLY relevant for output (.x<location> files)
// and ramp information (.rmps file);
// internally everything from 0..xmax!!!

double xmin;         
double tmin;        

// data and fields for on- or off-ramps

// choice_rmp=1

double x_rmp;
int    n_jumps_rmp;
double times_rmp[NBC], Q_rmp[NBC];
double flow_rmp;                 // assigned in "timestep" from Q_rmp[]


// choice_rmp>=2

int    n_rmps;                      // how many ramps
double x_rmps[NRMPMAX];
double dx_rmps[NRMPMAX];               // dx_rmp = input parameter, above
int    n_jumps_rmps[NRMPMAX];
double times_rmps[NRMPMAX][NBC];
double Q_rmps[NRMPMAX][NBC];
double flow_rmps[NRMPMAX];       // assigned in "timestep" from Q_rmps[][];


// choice_rmp==3: controlled ramps

int    n_rmps_c;
double x_rmps_c[NRMPMAX];
double dx_rmps_c[NRMPMAX];               // dx_rmp = input parameter, above
double s1_rmps_c[NRMPMAX];
double max_flow[NRMPMAX];
int    n_jumps_rmps_c[NRMPMAX];
double times_rmp_c[NBC], Q_rmp_c[NBC];
double times_rmps_c[NRMPMAX][NBC];
double Q_rmps_c[NRMPMAX][NBC];
double flow_rmps_c[NRMPMAX];       // assigned in "timestep" from Q_rmps[][];
double rmp_stor[NRMPMAX];
double Q_max;



// fields for time slices

int       n_slices;          // how many time slices
double time_slices[NBC];      // the actual times

// fields for x cuts

int   n_xcuts;             // how many x=const (loop) cross-sections
double x_xcuts[NBC];       // the actual locations

// fields for Dirichlet-BC

int n_jumps_l, n_jumps_r;
double times_l[NBC], rhoBC_l[NBC], QBC_l[NBC];
double times_r[NBC], rhoBC_r[NBC], QBC_r[NBC], vBC_r[NBC];


// ranges for debugging info if test_timestep is on

int    show          = false;
int it_show_min,it_show_max;
int ix_show_min,ix_show_max;
 


/*******************************************************************/
/* Density dependent fields and initial/boundary conditions
   calculated in "calc_tables" from the parameters above */
/*******************************************************************/

double Qmax;               // maximum equilibrium flow
double rhoQmax;            // rho value where euqil. flow=max

double Btab[NDV+1];           // universal dependence B(delta_V) of the 
                           // interaction term of the velocity eq.

double sqrtAtab[NRHO+1];   // A = Dimensionless measure for acceleration 
                           // fluctuations,definied by D(v)=2 A v^2/tau  
double Arhomax;            // A(rhomax)
double veqtab[NRHO+1];     // equilibrium velocity
double rho_congtab[NRHO+1];// Congested branch of the inverse Qe(rho) funct
double rho_freetab[NRHO+1];// Free branch of the inverse Qe(rho) funct
double rho_vtab[NRHO+1];   // inverse of the ve(rho) relation
double v0factab[NRHO+1];   // Reduction factor for desired velocity in the
                           // high-density correction or for KKK model
double high_denstab[NRHO+1];  
                        // in [0,1]; degree of the high-density correction 


// ***********************************************************************
// vars and fields for the t-dependence of the model parameters
// through varying truck fractions
// ***********************************************************************

class Trucks
{
 public:
  Trucks();
  int      varTruck;
  void     getTruckFrac(char namepar[]);  // defines fracTruckArr, timesArr
  void     getTruckPar(char namepar[]);   // defines v0Truck,TrTruck
  void     calc_corrections(double t);    // calculates v0TruckCorr,..
  double   v0TruckCorr;
  double   TrTruckCorr;
  double   rhomaxTruckCorr;

 private:
  double   v0Truck;
  double   TrTruck;
  double   rhomaxTruck;
  double   fracTruckArr [NBC+1];
  double   timesArr [NBC+1];
  int      nArr;    // field length of fracTruckArr, timesArr;
};

Trucks trucks=Trucks();

/***********************************************************************/
/*       Additional variables for macrosscopic Bando Model             */
/***********************************************************************/

class OVMparams
{
public:
  OVMparams();
  void getOVMparams (char namepar[]);
  double V0;// Parameters of the equilibrium velocity density relation 
  double V1;
  double c1;
  double c0;
};

OVMparams ovm = OVMparams();

/***********************************************************************/
/*       Additional variables for macroscopic IDM                      */
/***********************************************************************/

class IDMparams
{
public:
  IDMparams();
  void getIDMparams (char namepar[]);
  void calc_eq( double vwtab[],  double veqtab[]);
  double b;       // desired deceleration
  double s0;      // net distance maintained in standing traffic
  double delta;   // exponent
};

IDMparams idm = IDMparams();

/***********************************************************************/
/*       Additional variables for KKK and KKL models                   */
/***********************************************************************/


KKLparams kkl = KKLparams();

// initialize FPE* properly after reading input params

FPE fpe = FPE(); 
FPE_innov fpe_innov = FPE_innov();
VMM vmm = VMM();
SGM sgm = SGM();
Mac3phases mac3phases=Mac3phases();

/*******************************************************************/
/**************** Functions ****************************************/
/*******************************************************************/

void   add_rmp (double S1[], double S2[],  double rho[], double Q[],
              double x_rmp, double dx_rmp, double Q_rmp, int gauss);


void   advance_dt (int choice_method,
		double       u1[], double       u2[],
                const double F1[], const double F2[],
                const double S1[], const double S2[],
                int choice_BC,
                double    u1left,  double    u2left, 
                double    u1right, double    u2right,
                double dt, 
                double D1[], double D2[], 
                int show);

void   calc_eq_GKT(double rho, double vw, double A, double &veq);

void   boundary_vals(double field[], int nx, 
               int choice_BC, double left, double right);

void   calc_tables(char namepar[]);

void   calc_rhs (int choice_model, bool downwind_diff, 
              double rho[],double Q[],
              double F1[], double F2[],
              double S1[], double S2[],
              int show_calc_rhs);

// only test for isolated GKT calculation for exporting code

void   calc_rhsGKT (int choice_model, bool downwind_diff, 
              double rho[],double Q[],
              double F1[], double F2[],
              double S1[], double S2[],
              int show_calc_rhs);

double  cont_rmpflow(double flow_emp,double flow_sond1, double flow_crit, 
		     double dens_sond, double *stor_ptr);

void   control_test(char* func_name, int& flag);

void   cprhoQ (const double rho[], const double Q[],
              double rhocp[], double Qcp[]);

void   dyn_v0(double rho[]);

int    input(char namepar[] );

void   get_array (char* fname, int & n_array, 
          double array[]);
void   get_array (char* fname, int & n_array, 
          double array1[], double array2[]);
void   get_array (char* fname, int & n_array, 
          double array1[], double array2[], double array3[]);

void  get_array (char* fname, int & n_array, 
      double array1[], double array2[], double array3[], double array4[]);

void  get_dyn_v0(char namepar[], double x_v0[], double v0_x[],
		double delta_v0_x[],double thrsh_x[]);

void   get_xt_start(char namepar[], double& xmin, double& tmin);

void   getvar(FILE *fp, int *pint);
void   getvar(FILE *fp, double *pdouble);

void   get_xdependence(char namepar[], char varname[], 
                     double var_inp, double var_field[]);

void   get_T_xtdependence (char namepar[], double Tr_inp, 
                       double Tr_loc[], bool& flag_Txt);

void   initialize(double rho[], double Q[]);

void McCormack (  double       u1[], double       u2[],
                const double F1[], const double F2[],
                const double S1[], const double S2[],
                int choice_BC,
                double     u1left, double     u2left,
                double    u1right, double    u2right,
                double dt,
                double D1[], double D2[],
                int    show_McCormack);

void   shift_in_x (double antic_factor, 
                const double v[],  const double f[], double fshifted[]);


void shift_in_x_mca (const double rho_t[], 
                     const double f[], double fshifted[]);

void show_advance(const double u1start[], const double u2start[],
                  const double u1step[], const double u2step[], 
                  const double F1[], const double F2[], 
                  const double S1[], const double S2[]); 

void   timestep (int it, double rho[], double Q[], double a[]); 

void Upwind (   double       u1[], double       u2[],
                const double F1[], const double F2[],
                const double S1[], const double S2[],
                int choice_BC,
                double     u1left, double     u2left,
                double    u1right, double    u2right,
                double dt,
                double D1[], double D2[],
                int    show);

void  write_array (char* fname, int & n_array, 
      double array[]);

void   write_results(char npar[], int it, 
              const double rho[], const double Q[], const double a[]);

void   write_cross_section (int t_const, double t_sec, double x_meter, 
              const double rho[], const double Q[], const double a[]);

void simmaster ();

/**************************** end of declarations  **********************/

#endif //MASTER_H


