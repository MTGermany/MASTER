#ifndef GENERAL_H
#define GENERAL_H

#include <stdlib.h>

#include <string>
using namespace std;

/* General definitions for my C++ programs */
/* Martin Treiber */

#define PI acos(-1.)     
#define TRUE 1
#define FALSE 0
#define HUGE_VALUE 1.e6
#define TINY_VALUE 1.e-6

const int NVALS = 10000;

//warum mit static var???
//geht nicht auch a*a ????
#define SQR(x) ((x)*(x))
//static double sqrarg;
//#define SQR(a) ((sqrarg=(a))==0.0 ? 0.0 : sqrarg*sqrarg)

// Achtung: darf hier nix direkt definieren, sonst Fehler "multiple def"
// beim Linken (da general.h oft eingebunden)
// ausser: Explizit "inline" gegeben!

/*
inline double max(double x1, double x2){ return( (x1>x2) ? x1 : x2);}
inline double min(double x1, double x2){ return( (x1<x2) ? x1 : x2);}
inline int max (int x1, int x2){ return( (x1>x2) ? x1 : x2);} 
inline int min (int x1, int x2){ return( (x1<x2) ? x1 : x2);}
*/

double max(double x1, double x2);
double min(double x1, double x2);
int max (int x1, int x2);
int min (int x1, int x2);
//#######################################################
// arne, february 2006
// better parsing and convert methods using strings etc.

string itos(int i);
string dtos(double d);
double parseDouble(string str);
double parseDoubleWithSlash(string str);
bool parseDouble(string str, double &d);
int parseInt(string str);
unsigned int parseUnsignedInt(string str);
bool parseInt(string str, int &i);
bool parseBool(string str);
bool parseBool(string str, bool &b);

string trim_whitespace(string line);

bool fileExists(string filename);
//#######################################################
//double maxsmooth (double x1, double x2, double dx);  // in VW.cc
//double minsmooth (double x1, double x2, double dx);

int error (string errmsg);
double intp (const double tab[], int n, double x, double xmin, double xmax);
double intpextp (const double x_vals[], const double y_vals[], int n, double x);
double intpextp (const double x_vals[], const double y_vals[], double x,
                 int imin, int imax, bool reverse);
double intpextp (const double tab[], int n,double x, double xmin, double xmax);
void   filecheck(FILE *fp, const char filename[]);

double integrate(const double f[], int imin, int imax);
double normDist(const double w[], int imin, int imax, double xmin, double dx);
double avg(const double f[], int imin, int imax);
double avgDist(const double w[], int imin, int imax, double xmin, double dx);
double variance(const double f[], int imin, int imax);
double varianceDist(const double w[], int imin, int imax, double xmin, double dx);
double skewness(const double f[], int imin, int imax);
double skewnessDist(const double w[], int imin, int imax, double xmin, double dx);
double kurtosis(const double f[], int imin, int imax);
double kurtosisDist(const double w[], int imin, int imax, double xmin, double dx);
double scalarProd(const double f[], const double g[], int imin, int imax);
double corr(const double f[], const double g[], int imin, int imax);
void shift (const double f[], int n, int ishift, double fshifted[],int& imin, int& imax);


//double roundDouble(double doValue, int nPrecision);

/* Table of the error function from x=-3...3 */

const double XERFMIN  = -3.;
const double XERFMAX  = 3.;
const int NERFTAB  = 12;          // NERFTAB+1 = Number of points in the
                                  // table "erftab" for the error function

const double erftab[] = {0.0013, 0.0062, 0.0228, 0.0668, 0.1587, 0.3085,
    0.5, 0.6915, 0.8413, 0.9332, 0.9772, 0.9938, 0.9987};

const double sqrt2 = 1.414213562;
#endif // GENERAL_H
