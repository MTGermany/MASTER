
/*******************************************************************/
/*  General selfmade c and c++ functions 
    general io routines (get_array*, write_array* ...) in io.cc */
/*******************************************************************/
//C
#include <stdio.h>
#include <math.h>

//C++
#include <iostream>
#include <fstream>
#include <sstream>

#include "general.h"





//#####################################################################
//  function: convert integer to string 
//#####################################################################

string itos(int i)
{
  stringstream s;
  s << i;
  return s.str();
}


//#####################################################################
//  function: convert double to string   
//#####################################################################


string dtos(double d)
{
  stringstream s;
  s << d;
  return s.str();
}



//#####################################################################
//  function: parse double from string
//#####################################################################

double parseDouble(string str)
{
  double d;
  istringstream iss(str);
  iss >> d;
  return d;
}

double parseDoubleWithSlash(string str){
	string::size_type loc = str.find("/", 0 );
	if( loc == string::npos ) return parseDouble(str);
	
	string p = str.substr(0, loc);
	string q = str.substr(loc+1, str.length());
	//cout<<"p = "<<p<<", q = "<<q<<endl;
	double d1;
	istringstream iss1(p);
	iss1 >> d1;
	double d2;
	istringstream iss2(q);
	iss2 >> d2;
	return d1/d2;
}
//NOCH TESTEN: char* strings
bool parseDouble(string str, double &d)
{
  d=0;
  istringstream iss(str);
  iss >> d;
  //error treatment: 
  char ch;
  if ((!iss) || (iss >> ch)) return(false);
  return(true);
}

//#####################################################################
//  function: parse double from string
//#####################################################################

int parseInt(string str)
{
  int i;
  istringstream iss(str);
  iss >> i;
  return i;
}

unsigned int parseUnsignedInt(string str)
{
  unsigned int i;
  istringstream iss(str);
  iss >> i;
  return i;
}

//noch testen!!!
bool parseInt(string str, int &i)
{
  i=0;
  istringstream iss(str);
  iss >> i;
  //error treatment: 
  char ch;
  if((!iss)||( iss >> ch )) return(false);
  return(true);
}

//#####################################################################
//  function: parse double from string
//#####################################################################

bool parseBool(string str)
{
   return (str=="true" || str=="True" || str=="TRUE");
}

bool parseBool(string str, bool &b)
{
  b=false;
  if(str=="true" || str=="TRUE" || str=="True")
    {
      b=true;
      return true;
    }
  if(str=="false" || str=="FALSE" || str=="False")
    {
      b=false;
      return true;
    }
  return false;
}

string trim_whitespace(string line){
  if(line.empty()) return line;
  unsigned int i;
  for(i=0; i<line.size(); i++){
    if(line[i] != ' ' && line[i] != '\t' ) break;
  }
  return line.substr(i);
}

//#####################################################################
//  function: check file exists?
//#####################################################################

bool fileExists(string filename)
{
  ifstream infile(filename.c_str(), ios::in);
  if(!infile)return false;
	else{
		infile.close();
		return true;
	}
}

//#####################################################################
//  function  intp   
//#####################################################################

double intp (const double tab[], int n,double x, 
	     double xmin, double xmax)
   /* intp interpolates the array tab with n 
 (before aug 04: n+1!!!) equidistant points
      in [xmin,xmax] at the location x; an error message is produced on
      attempting extrapolation */
{
  int nloc=n-1; //!!
  double intp_value;
  double ir   = nloc*(x-xmin)/(xmax-xmin);
  int    i    = (int) ir;
  double rest = ir-i;
  if ((i>=0) && (i<=nloc-1))  intp_value =  (1-rest) * tab[i] + rest*tab[i+1];
  else if (i==nloc) intp_value = tab[nloc];
  else {
    cout << "intp: index i = "<<i<<" (ir="<<ir<<") out of range\n";
    exit(1);
  }

  return(intp_value);
}

//#####################################################################
//  function  intpextp   
//#####################################################################

double intpextp (const double tab[], int n,double x, 
		 double xmin, double xmax)
   /* As intp, but extrapolates the constant
      boundary values instead of an error message */
{
  int nloc=n-1; //!!
  double intp_value;
  double ir   = nloc*(x-xmin)/(xmax-xmin);
  int    i    = (int) ir;
  double rest = ir-i;
  if ((i>=0) && (i<nloc-1))  intp_value =  (1-rest) * tab[i] + rest*tab[i+1];
  else if (i>=nloc-1) intp_value = tab[nloc];
  else intp_value = tab[1];

  return(intp_value);
}


double intpextp (const double x_vals[], const double y_vals[], 
                        int n, double x)
   /* Inter- and extrapolation of tabulated
      functions given by an array  x_vals with INCREASING 
      (but nonequidistant) values, 
      and the associated y values y_vals.
      x is the value for which inter(extra-)polation is sought.
      On extrapolation, the left value y_vals[0]
      (right value y_vals[n]) is outputted if x<x_vals[0] (x>x_vals[n]).
       NOT TIME OPTIMIZED  
   */
{
    int i=0;
    int nloc=n-1; //!!
    double intp_value;

    while( (x_vals[i] <= x) && (i<nloc)) i++;

    if(i==0) intp_value = y_vals[0];
    else if ( (i==nloc) && (x > x_vals[i])) 
             intp_value = y_vals[nloc];
    else if (fabs(x_vals[i]-x_vals[i-1])<TINY_VALUE)
             intp_value = y_vals[i];
    else     intp_value = y_vals[i-1] + (y_vals[i]-y_vals[i-1]) 
               * (x - x_vals[i-1])/(x_vals[i]-x_vals[i-1]);

    return(intp_value);
}

double intpextp (const double x_vals[], const double y_vals[], double x,
    int imin, int imax, bool reverse)
   /* Version using only the index range imin ... imax and extrapolating
      constant values otherwise; reverse=true means that the array x_vals
      has x values in decreasing order. NOT TIME OPTIMIZED
 */
{
    double intp_value;

    int i=imin;
    if(reverse)
      while( (x_vals[i] >= x) && (i<imax)) i++;
    else
      while( (x_vals[i] <= x) && (i<imax)) i++;

    if(i==imin)                                      // left extrapolation
              intp_value = y_vals[imin];            
    else if (i==imax)                                // right extrapolation
             intp_value = y_vals[imax];
    else if (fabs(x_vals[i]-x_vals[i-1])<TINY_VALUE) // same x values
             intp_value = y_vals[i];
    else                                             // interpolation
            intp_value = y_vals[i-1] + (y_vals[i]-y_vals[i-1])  
              * (x - x_vals[i-1])/(x_vals[i]-x_vals[i-1]);
    return(intp_value);
}
/********************************************************************/
//  arne (2005) 
/********************************************************************/

double max(double x1, double x2){ return( (x1>x2) ? x1 : x2);}
double min(double x1, double x2){ return( (x1<x2) ? x1 : x2);}
int max (int x1, int x2){ return( (x1>x2) ? x1 : x2);} 
int min (int x1, int x2){ return( (x1<x2) ? x1 : x2);}



/********************************************************************/
/*  RoundDouble arne (2005)                                          */
/********************************************************************/
/*
double roundDouble(double doValue, int nPrecision)
{
  //this version found on web:
  static const double doBase = 10.0;
  double doComplete5, doComplete5i;
  
  doComplete5 = doValue * pow(doBase, static_cast<double>(nPrecision + 1));
  
  if(doValue < 0.0)
    doComplete5 -= 5.0;
    else
      doComplete5 += 5.0;
  
  doComplete5 /= doBase;
  modf(doComplete5, &doComplete5i);
  return ( doComplete5i / pow(doBase, static_cast<double>(nPrecision)) );
}
*/

/********************************************************************/
/*  function  error                                                 */
/********************************************************************/

int error (string errmsg)
{ 
  fprintf(stderr,"%s\n",errmsg.c_str()); 
  exit(1);
}


/********************************************************************/
/*  function  filecheck                                             */
/********************************************************************/

void filecheck(FILE *fp, const char filename[])
{
  char message[80];
  if(fp==NULL)
  {
    sprintf(message,"fopen: Error while opening the file %s\n",filename);
    error(message);
  }
}



//########################################################################
//arne: habe von double zu double umgetypt
double integrate(const double f[], int imin, int imax)
//########################################################################
    // trapezoid rule; equidistant points with distane 1 assumed
  {
    int i;
    double result=(f[imin] + f[imax]) / 2.;
    for (i=imin+1; i<=imax-1; i++) result += f[i];
    return(result);
  }





//########################################################################
  double normDist(const double w[], int imin, int imax, double xmin, double dx)
//########################################################################

    // Normalize or test normalization of Distribution function w[]
    // with nonzero values from imin to imax, 
    // and abszissa values x[i]=xmin + (i-imin)*dx
  {
    double result=0.;
    for(int i=imin; i<=imax; i++)
    {
      result += w[i]*dx;
    }
    return(result);
  }

//########################################################################
  double avg(const double f[], int imin, int imax)
//########################################################################

    // expectation value of field of values
  {
    const int n = imax-imin+1; 
       if(n<1) error("avg: field has zero length!");
    int       i;
    double result=0.;
    for(i=imin; i<=imax; i++) result+=f[i];
    result/=n;
    return(result);
  }
 
//########################################################################
  double avgDist(const double w[], int imin, int imax, double xmin, double dx)
//########################################################################

    // expectation value of normalized Distribution function w[]
    // with nonzero values from imin to imax, 
    // and abszissa values x[i]=xmin + (i-imin)*dx
  {
    double result=0.;
    for(int i=imin; i<=imax; i++)
    {
      double x = xmin+(i-imin)*dx;
      result += x * w[i]*dx;
    }
    return(result);
  }

//########################################################################
  double variance(const double f[], int imin, int imax)
//########################################################################
  {
    const int n  = imax-imin+1; 
       if(n<2) error("var: field has length < 2!");
    int   i;
    double mu   = avg(f,imin,imax);
    double sqrsum = 0.;
    for(i=imin; i<=imax; i++) sqrsum += SQR(f[i]);

    double result = ( sqrsum - n*SQR(mu)) / (n-1); // emp. variance!!
    return(result);
  }

//########################################################################
  double varianceDist(const double w[], int imin, int imax, double xmin, double dx)
//########################################################################

    // expectation value of normalized Distribution function w[]
  {
    double result=0.;
    double mu = avgDist(w,imin,imax,xmin,dx);
    for(int i=imin; i<=imax; i++)
    {
      double x = xmin+(i-imin)*dx;
      result += SQR(x-mu) * w[i]*dx;
    }
    return(result);
  }

//########################################################################
  double skewness(const double f[], int imin, int imax)
//########################################################################
  {
    const int n  = imax-imin+1; 
       if(n<2) error("var: field has length < 2!");
    int   i;
    double mu   = avg(f,imin,imax);
    double sigma   = sqrt( variance(f,imin,imax));
    double erwpow3 = 0.;
    for(i=imin; i<=imax; i++) erwpow3 += pow( (f[i]-mu), 3)/n;

    double result = erwpow3/pow(sigma,3);
    return(result);
  }


//########################################################################
  double skewnessDist(const double w[], int imin, int imax, double xmin, double dx)
//########################################################################

    // expectation value of normalized Distribution function w[]
  {
    double result=0.;
    double mu = avgDist(w,imin,imax,xmin,dx);
    double sigma = sqrt( varianceDist(w,imin,imax,xmin,dx) );
    for(int i=imin; i<=imax; i++)
    {
      double x = xmin+(i-imin)*dx;
      result += pow( (x-mu),3) * w[i]*dx;
    }
    result/=pow(sigma,3);
    return(result);
  }


//########################################################################
double kurtosis(const double f[], int imin, int imax)
//########################################################################
  {
    const int n  = imax-imin+1; 
       if(n<2) error("var: field has length < 2!");
    int   i;
    double mu   = avg(f,imin,imax);
    double sigma   = sqrt( variance(f,imin,imax));
    double erwpow4 = 0.;
    for(i=imin; i<=imax; i++) erwpow4 += pow( (f[i]-mu), 4)/n;

    double result = erwpow4/pow(sigma,4) - 3.;
    return(result);
  }



//########################################################################
  double kurtosisDist(const double w[], int imin, int imax, double xmin, double dx)
//########################################################################

    // expectation value of normalized Distribution function w[]
  {
    double result=0.;
    double mu = avgDist(w,imin,imax,xmin,dx);
    double sigma = sqrt( varianceDist(w,imin,imax,xmin,dx) );
    for(int i=imin; i<=imax; i++)
    {
      double x = xmin+(i-imin)*dx;
      result += pow( (x-mu),4) * w[i]*dx;
    }
    result = result / pow(sigma,4) - 3;
    return(result);
  }



//########################################################################
  double scalarProd(const double f[], const double g[], int imin, int imax)
//########################################################################
  {
    int i;
    double result=(f[imin]*g[imin] + f[imax]*g[imax]) / 2.;
    for (i=imin+1; i<=imax-1; i++) result += f[i]*g[i];
    return(result);
  }



//########################################################################
  double corr(const double f[], const double g[], int imin, int imax)
//########################################################################

  // general correlation of two fields going from 0 to n
  // corr_f(dt) = <ferw0_t, ferw0_{t+dt}> / sqrt(||ferw0_t|| ||ferw0_{t+dt}||)
  // with the norm ||f||^2 = <f,f>, ferw0 = f - <f>
  // and <f,g> = integral from f*g in considered t range

  // integrate and scalarProd do the dirty work: 
  // integration with the trapezoid rule
  {
    if(imax<=1) error("corr: correlation of function with <=2 points attempted");
    if((imin<0) || (imax >NVALS)  ) 
      error("corr: array out of bounds (imax>NVALS) or imin<0)");
    if(imin>=imax) error("corr: imin>=imax");

    int       i;
    const int n=imax-imin+1;

    //double     ferw0[NVALS+1];
    //   double     gerw0[NVALS+1];
    double     *ferw0; ferw0 = new double[imax+1];
    double     *gerw0; gerw0 = new double[imax+1];

    double     erwf  = integrate(f,imin, imax) / n;
    double     erwg  = integrate(g,imin, imax) / n;

    for (i=imin; i<=imax; i++) ferw0[i] = f[i]-erwf;
    for (i=imin; i<=imax; i++) gerw0[i] = g[i]-erwg;

    double     normf = sqrt(scalarProd(ferw0,ferw0,imin, imax));
    double     normg = sqrt(scalarProd(gerw0,gerw0,imin, imax));

    double     denom = normf*normg;
      if(denom==0) error("corr: At least one of the functions is zero");

    double     result = scalarProd(ferw0,gerw0,imin, imax)/denom;

    delete[] ferw0;
    delete[] gerw0;
    return(result);
  }


//########################################################################
  void shift (const double f[], int n, int ishift, double fshifted[],  
              int& imin, int& imax)
//########################################################################
  {
   
    // f[] has (n+1) entries from 0 to n !!!
    // min/max position relative  to unshifted field

    int i;
    imin = (ishift>0) ? 0 : -ishift;
    imax = (ishift>0) ? n-ishift : n;
    for (i=imin; i<=imax; i++) fshifted[i] = f[i+ishift];

    /*    cout << "after shift: n="<<n<<"f[0]= "<<f[0]
              <<" f[n]= "<<f[n]  
              <<" imin= "<<imin  
              <<" imax= "<<imax  
	 << endl;
    */
   }
