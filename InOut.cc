

//(jun 04) ACHTUNG: InOut gibt nun n_array als Zahl der Eintraege raus!
// max. Index also n_array-1  !!

//arne (may 05): habe <strstream.h> (VERALTET) ausgestauscht gegen sstream
//die klasse heisst nun istringstream anstatt istrstream

// c++ 
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <cstring>
#include <cstdlib>

 
using namespace std;

#include "InOut.h"


InOut::InOut(){;}


bool InOut::fileExists(const char* fname){
  ifstream  infile (fname, ios::in);

  // formulate so cumbersom since nature of infile
  // (bool? file descriptor?) obscure

  if(!infile){return false;}
  else return true;
}



int InOut::getNumberOfLines(const char* fname)
{
  cout<<"in InOut.getNumberOfLines:: " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if(!infile)
    {
      cerr << "Error opening file " << fname << " for reading" << endl;
      exit(-1);
    }
  
  int i=0;
  char line[LINEMAX];
  //while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  //achtung darf nicht aufhoeren bei leerzeile...
  while ( (infile.getline(line,LINEMAX))  )
    {
      if( is_not_white(line,strlen(line)) && is_data_line(line) )
	{
	  i++;
	}
    }
  infile.close();
  cout << "number of lines = " <<i<<endl;
  return(i);
}

int InOut::getNumberOfCols(const char* fname)
{
  cout<<"in InOut.getNumberOfCols:: " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if(!infile)
    {
      cerr << "Error opening file " << fname << " for reading" << endl;
      exit(-1);
    }
  
  char line[LINEMAX];
  bool dataLineReached=false;

  //while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  //achtung darf nicht aufhoeren bei leerzeile...
  while ( (infile.getline(line,LINEMAX))&&(!dataLineReached)  ){
    if( is_not_white(line,strlen(line)) && is_data_line(line) ){
      dataLineReached=true;
    }
  }
  infile.close();
  int ncol=0;
  if(!dataLineReached){
    cerr<<"no data line found in "<<fname<<" return zero columns"<<endl;
    return(0);
  }
  else{
    istringstream linestream (line);
    double proforma;
    while(linestream){
      linestream >> proforma;
      ncol++;
    }
    ncol--;
    cout << "number of cols = " <<ncol<<endl;
  }
  return(ncol);
}





void  InOut::get_array(char* fname, int & n_array, double array1[])
{
  cout << "in InOut.get_array (double array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  //while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  while ( infile.getline(line,LINEMAX) )
  {
    if( is_data_line(line) && is_not_white(line,strlen(line)) ){
      prepare_line(line);
      istringstream linestream (line);
      linestream >> array1[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}

void  InOut::get_array (char* fname, int & n_array, int array1[]){
  cout << "in InOut.get_array (int array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      linestream >> array1[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}

void  InOut::get_col (const char* fname, int col,  int & n_array, int array1[]){
  cout << "in InOut.get_col(int array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  //while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  while ( (infile.getline(line,LINEMAX))) 
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<col; icol++){linestream >> proforma;}
      linestream >> array1[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}

void  InOut::get_col (const char* fname, int col,  int &n_array, double array1[]){
  cout << "in get_col(double array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  //  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  
  while ( (infile.getline(line,LINEMAX))) 
    {
		
      if(is_data_line(line)){
	prepare_line(line);
	//if(i<4){cout<<"InOut.line="<<line<<endl;}
	istringstream linestream (line);
	for (int icol=1; icol<col; icol++){linestream >> proforma;}
	//	char ch;
	linestream>>array1[i];
	//error treatment: arne feb 08
	//testweise !!!
	if (!linestream){
	  cerr<<"error while parsing for double ... i="<<i<<" --> "
              <<array1[i]<<" ... set to 0 ..."<<endl;
          array1[i]=0;
	} 
	//cout<<"i="<<i<<" InOut.array1[i]="<<array1[i]<<" n_array="<<n_array<<endl;
	i++;
      }
    }
	
  infile.close();
  n_array = i;
  cout << "InOut.get_col: col="<<col<<" n_array = " << n_array << endl;
}

void  InOut::getChar_col (const char* fname, int col, int & n_array, char array1[]){
  cout << "in InOut.getChar_col( char array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<=col; icol++){linestream >> proforma;}
      array1[i]=proforma[0];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}



//arne 9-7-07
void  InOut::get_col (const char* fname, int col,  int & n_array, string array1[]){
  cout << "in InOut.get_col(int array): " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[LINEMAX];
  char proforma[LINEMAX]; 
  while ( (infile.getline(line,LINEMAX))) 
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      for (int icol=1; icol<col; icol++){linestream >> proforma;}
      linestream >> array1[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}




void  InOut::get_array (char* fname, int & n_array, 
      double array1[], double array2[])
{
  cout << "in InOut.get_array (double, double): " 
       << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;

  //  while (infile >> array1[i] >> array2[i]) i++;

  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
      linestream >> array1[i] >> array2[i];
      cout << "array1[" << i <<"] = " << array1[i] 
             << " array2[" << i <<"] = " << array2[i] << endl;
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}


void  InOut::get_array (char* fname, int & n_array, 
      double array1[], double array2[], double array3[])
{
  cout << "in InOut.get_array(double,double,double): " 
       << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;

  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
        istringstream linestream (line);
        linestream >> array1[i] >> array2[i] >> array3[i];
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout << "n_array = " << n_array << endl;
}

void  InOut::get_array (char* fname, int & n_array, 
      int array1[], double array2[], double array3[])
{
  cout << "in InOut.get_array (int, double, double): " 
       << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;

  char line[LINEMAX];
  while ( (infile.getline(line,LINEMAX))  && (is_not_white(line,strlen(line))) )
  {
    if(is_data_line(line)){
      prepare_line(line);
      istringstream linestream (line);
       linestream >> array1[i] >> array2[i] >> array3[i];
      cout <<"InOut::get_array(*,*,int[], ...):"
	   <<" line="<<line<<endl
	   <<" i="<<i<<" array1[i]="<<array1[i]<<endl;
      i++;
    }
  }
  infile.close();
  n_array = i;
  cout <<"n_array = " << n_array << endl;
}




//################################################################
//  function  get_array2d 
//################################################################

void  InOut::get_array2d (char* fname, 
                  double unitTimeInDat_s, double unitSpaceInDat_m,
                  int & nt, int & nx, double & dt, double & dx, 
                  double array1[NXMAX][NTMAX+1],
                  double array2[NXMAX][NTMAX+1],
                  double array3[NXMAX][NTMAX+1])

{
  cout << "in get_array2d: " << "fname = "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  const int MAXSTR=500;
  char line[MAXSTR];
  double xvals[NXMAX];
  double tvals[NTMAX];

  int i=0;
  while (!(infile.eof())){

    int j=0;
    if (i>=NTMAX){
      cerr << "InOut.get_array2d: error: time index (second coordinate) i="<<i
	   <<" >=  NTMAX="<<NTMAX<<endl;
      exit(-1);
    }
    while (   (infile.getline(line,MAXSTR))  
           && (is_not_white(line,strlen(line))) 
          )

    {
      if( (line[0] != COMMENTCHAR) && 
          (line[0] != COMMENTCHAR2) &&
          (strlen(line)>0) )
      {
	//	cout <<"i="<<i<<" j="<<j<<endl;

        std::istringstream linestream (line);
        linestream >>  xvals[j] >>  tvals[i] 
                   >> array1[j][i] >> array2[j][i] >> array3[j][i];
        if(false)cout << "array1[" << i <<"][" << j <<"] = " << array1[j][i] << endl;
        j++;
      }
      nx = j-1; //!!
      nx = j;
    }
    i++;
  }

  infile.close();
  nt = i-2; //!!
  nt = i-1;
  dx = (nx>=2) ? unitSpaceInDat_m*(xvals[nx-1] - xvals[nx-2]) : xvals[nx-1];
  dt = (nt>=2) ? unitTimeInDat_s*(tvals[nt-1] - tvals[nt-2]) : xvals[nt-1];

  //  for (i=0; i<=nt; i++) cout << "tvals[i]="<<tvals[i]<<endl;
  cout << "nt = " << nt << "; nx = " << nx << endl;
  cout << "dt = " << dt << "; dx = " << dx << endl;
  cout << "array1[" << nx-1 <<"][" << nt-1 <<"] = " << array1[nx-1][nt-1] << endl
       << "array1[" << 0 <<"][" << 0 <<"] = " << array1[0][0] << endl
       << "array2[" << nx-1 <<"][" << nt-1 <<"] = " << array2[nx-1][nt-1] << endl
       << "array3[" << nx-1 <<"][" << nt-1 <<"] = " << array3[nx-1][nt-1] << endl;

  //exit(0);

}




void  InOut::write_array (char* fname, int nData, const double data_col1[],
			  char* name_col1){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", name_col1);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\n", data_col1[i]);
  }
  fclose(outfile);
} // InOut::write_array


void  InOut::write_array (char* fname, int nData, const double data_col1[],
const double data_col2[], char* name_col1, char* name_col2)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\n", name_col1, name_col2);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\n", data_col1[i], data_col2[i]);
  }
  fclose(outfile);
} // InOut::write_array


void  InOut::write_array (char* fname, int nData, const int data_col1[],
const int data_col2[], char* name_col1, char* name_col2)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\n", name_col1, name_col2);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%i\t%i\n", data_col1[i], data_col2[i]);
  }
  fclose(outfile);
} // InOut::write_array

void  InOut::write_array (char* fname, int nData, const double data_col1[],
const double data_col2[], const double data_col3[], 
char* name_col1, char* name_col2, char* name_col3)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\n", name_col1, name_col2, name_col3);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\t%f\n", data_col1[i], data_col2[i], data_col3[i]);
  }
  fclose(outfile);

} // InOut::write_array


//<new aug16>

void InOut::write_array (char* fname, int nData, const double data_col1[],
			 const double data_col2[], const double data_col3[],
			 char* titleString){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%1.7f\t%1.7f\n", 
	    data_col1[i], data_col2[i], data_col3[i]);
  }
  fclose(outfile);
} 

//</new>



void  InOut::write_array (char* fname, int nData, const double data_col1[],
const double data_col2[], const double data_col3[], const double data_col4[], 
char* name_col1, char* name_col2, char* name_col3, char* name_col4)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\n", name_col1, name_col2, 
        name_col3, name_col4);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%6.8f\t%f\t%6.8f\n", data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i]);
  }
  fclose(outfile);

} // InOut::write_array

void  InOut::write_array (char* fname, int nData, const double data_col1[],
const int data_col2[], const double data_col3[], const double data_col4[], 
char* name_col1, char* name_col2, char* name_col3, char* name_col4)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\n", name_col1, name_col2, 
        name_col3, name_col4);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%i\t%f\t%6.8f\n", data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i]);
  }
  fclose(outfile);
} // InOut::write_array


//<new nov17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i]);
  }
  fclose(outfile);
} 
//</new>


void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\t%f\t%f\t%f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i]);
  }
  fclose(outfile);

} // InOut::write_array

void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const int data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%i\t%f\t%f\t%f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i]);
  }
  fclose(outfile);

} // InOut::write_array




//<new nov17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],
		     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
      //fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i]);
  }
  fclose(outfile);
} 
//</new>



void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%f\t%f\t%f\t%f\t%f\t%f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i]);
  }
  fclose(outfile);
} 


//<new nov17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
		     const double data_col5[],const double data_col6[],
		     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i]);
  }
  fclose(outfile);
} 
//</new>




void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i]);
  }
  fclose(outfile);

} // InOut::write_array


//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],
                     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i]);
  }
  fclose(outfile);

} // InOut::write_array





void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7, name_col8);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i]);
  }
  fclose(outfile);

} // InOut::write_array

//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],const double data_col8[],
                     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i]);
  }
  fclose(outfile);

} // InOut::write_array







void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8, 
                     char* name_col9)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7, name_col8, 
        name_col9);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i]);
  }
  fclose(outfile);

} // InOut::write_array

//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],const double data_col8[],
                     const double data_col9[],
			  char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i]);
  }
  fclose(outfile);

} // InOut::write_array









void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8, 
                     char* name_col9, char* name_col10)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
        name_col1, name_col2, 
        name_col3, name_col4, 
        name_col5, name_col6, 
        name_col7, name_col8, 
        name_col9, name_col10);
  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i], data_col10[i]);
  }
  fclose(outfile);

} // InOut::write_array

//<new mar17>
void  InOut::write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
		     const double data_col7[],const double data_col8[],
                     const double data_col9[],const double data_col10[],
		     char* titleString)
{
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "#%s\n", titleString);

  for(int i=0; i<nData; i++){
    fprintf(outfile, "%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
       data_col1[i], data_col2[i], 
       data_col3[i], data_col4[i], 
       data_col5[i], data_col6[i], 
       data_col7[i], data_col8[i], 
       data_col9[i], data_col10[i]);
  }
  fclose(outfile);

} // InOut::write_array


//##########################################################
// general-purpose function  write_array2d 
// nata  => nLines
// nCols => generalizes above versions of write_array
//##########################################################

void  InOut::write_array_gen (char* fname, int nLines, int nCols,  
                     const double data[][LINEMAX], const char* headers){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "%s\n", headers);
  for(int i=0; i<nLines; i++){
    char line[MAXSTR];
    sprintf(line,""); // even if compiler complains: This line *is* necessary; otherwise overflow!
    for(int icol=0; icol<nCols; icol++){
      //cout <<"test: i="<<i<<" icol="<<icol<<" data[icol][i]="<<data[icol][i]<<endl;
      sprintf(line,"%s%5.4f\t", line, data[icol][i]);
    }
    //cout<<"test: line="<<line<<endl;
    fprintf(outfile, "%s\n", line);
    
  }
}

void  InOut::write_array_gen (char* fname, int nLines, int nCols,  
                     const int data[][LINEMAX], const char* headers){
  FILE  *outfile;                  
  outfile = fopen(fname,"w");
  if(!outfile){cerr<<"InOut.write_array: Error: file "<<fname
		   <<" cannot be created or is not writable"<<endl;
    exit(-1);
  }

  fprintf(outfile, "%s\n", headers);
  for(int i=0; i<nLines; i++){
    char line[MAXSTR];
    sprintf(line,""); // even if compiler complains: This line *is* necessary; otherwise overflow!
    for(int icol=0; icol<nCols; icol++){
      sprintf(line,"%s%i\t", line, data[icol][i]);
    }
    fprintf(outfile, "%s\n", line);
    
  }
}



//########################################################################
//  function  get_array2d 
//########################################################################

/// reads in array[i][j] representing z(x_i,y_j) with
/// regular spacings 
/// x_i=xmin+i*dx, i=0 .. nx-1, 
/// y_j=ymin+j*dy, y=0 .. ny-1
/// from column col (>=3) of file fname (column count starts at 1); 
/// index i (first col) is the same in each block, 
/// index j (second col) varies every line; empty line between 
/// two blocks 
/// Output: 
/// xmin,dx,nx (calculated from first col)
/// ymin,dy,ny (calculated from second col)
/// array[i][j]

void  InOut::get_array2d_col (const char* fname, int col, 
		       double& xmin, double& dx, int& nx, 
		       double& ymin, double& dy, int& ny, 
		       double array[NXMAX][NYMAX]){
  cout << "in InOut.get_array2d: " 
       << "fname = "<<fname<<" col="<<col <<endl;
  if(col<3){
    cerr<<" InOut.get_array2d: Error: column index must be >=3"
	<<endl;
    exit(-1);
  }
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "InOut.get_array2d: Error opening file " << fname 
         << " for reading" << endl;
    exit(-1);
  }
  
  char line[LINEMAX];
  char proforma[LINEMAX]; 

  double xvals[NXMAX];
  double yvals[NYMAX];

  int i=0;
  int i_in_loop=0;
  while (!(infile.eof())){

    int j=0;
    while (   (infile.getline(line,LINEMAX))  
           && (is_not_white(line,strlen(line))) 
          ){
    if(is_data_line(line)){
      prepare_line(line);
        istringstream linestream (line);
        linestream >>  xvals[i] >>  yvals[j];
        for (int icol=3; icol<col; icol++){linestream >> proforma;}
        linestream >> array[i][j];
        cout << "array[" << i <<"][" << j <<"] = " << array[i][j] << endl;
        j++;
	i_in_loop=i;
      }
      ny = j;
    }
    i++;
    cout <<"after i++: i="<<i<<endl;
  }
  //nx = i-2; ???
  nx = i_in_loop+1;

  infile.close();

  xmin=xvals[0];
  ymin=yvals[0];
  dx = (nx>1) ? xvals[nx-1] - xvals[nx-2] : xvals[nx-1];
  dy = (ny>1) ? yvals[ny-1] - yvals[ny-2] : yvals[ny-1];

  if(true){
    cout <<"InOut.get_array2d_col successfully:\n";
    cout << "xmin = " << xmin<< "\tdx = " << dx << "\tnx = " << nx << endl;
    cout << "ymin = " << ymin<< "\tdy = " << dy << "\tny = " << ny << endl;
    cout << "array[" << 0  <<"][" << 0  <<"] = " << array[0][0]    << endl
	 << "array[" << 0  <<"][" << ny-1 <<"] = " << array[0][ny-1]  << endl
	 << "array[" << nx-1 <<"][" << 0  <<"] = " << array[nx-1][0]  << endl
	 << "array[" <<nx-1<<"]["<< ny-1 <<"] = " << array[nx-1][ny-1] << endl;
    }

}

//##########################################################
// general-purpose function  write_array2d 
// writes gnu-plottable file with col1=x value, col2=y value, col3=array=z(z,y) 
//##########################################################

// definition of arrays for calling write_array2d: see HINT in .h file
// see also write_array2d_fromarray1d

void  InOut::write_array2d (char* fname,
		           double xmin, double xmax, int nx,
		           double ymin, double ymax, int ny,
			    double array[][NYMAX],
			    char* titleString){

  if((nx>NXMAX)||(ny>NYMAX)){
    cerr <<"InOut::write_array2d:" 
                 <<" nx="<<nx<<" >= NXMAX="<<NXMAX<<" and/or "
	 <<" ny="<<ny<<"  >= NYMAX="<<NYMAX
	 <<" => Error!" << endl;
    exit(-1);
  }

  double dx=(xmax-xmin)/((double)(nx-1));
  double dy=(ymax-ymin)/((double)(ny-1));
  cout << "in write_array2d: "
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<" xmax="<<xmax <<" nx="<<nx<<" dx="<<dx
       <<" ymin="<<ymin<<" ymax="<<ymax <<" ny="<<ny<<" dy="<<dy <<endl;

  FILE  *outfile;                  

  outfile = fopen(fname,"w");

  fprintf(outfile, "#%s \n",titleString);

  for (int ix=0; ix<nx; ix++){
    for (int iy=0; iy<ny; iy++){
      double x=xmin+ix*dx;
      double y=ymin+iy*dy;
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\n", x,y,array[ix][iy]);
    }
    fprintf(outfile, "\n");        // newline mean new "sweep" for gnuplot
  }

  fclose(outfile);

} // end general-purpose write_array2d

//##########################################################
// general-purpose function  write_array2d_fromarray1d 
// writes gnu-plottable file with col1=x value, col2=y value, 
// col3=array=z(z,y) where z[i][j]=array1d[nx*i+j]
//##########################################################

// definition of arrays for calling write_array2d: see HINT in .h file

void  InOut::write_array2d_fromarray1d (char* fname,
		           double xmin, double xmax, int nx,
		           double ymin, double ymax, int ny,
			    double array1d[],
			    char* titleString){


  double dx=(xmax-xmin)/((double)(nx-1));
  double dy=(ymax-ymin)/((double)(ny-1));
  cout << "in write_array2d_fromarray1d: "
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<" xmax="<<xmax <<" nx="<<nx<<" dx="<<dx
       <<" ymin="<<ymin<<" ymax="<<ymax <<" ny="<<ny<<" dy="<<dy <<endl;

  FILE  *outfile;                  

  outfile = fopen(fname,"w");

  fprintf(outfile, "#%s \n",titleString);

  for (int ix=0; ix<nx; ix++){
    for (int iy=0; iy<ny; iy++){
      double x=xmin+ix*dx;
      double y=ymin+iy*dy;
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\n", x,y,array1d[nx*ix+iy]);
    }
    fprintf(outfile, "\n");        // newline mean new "sweep" for gnuplot
  }

  fclose(outfile);

} // end general-purpose write_array2d_fromarray1d

//#########################################################
//  special-purpose function  write_array2d (general purpose: above)
//#########################################################

  // definition of arrays for calling write_array2d: see HINT in .h file

void  InOut::write_array2d (char* fname, double dxout, double dtout, 
			    double xmin, double xmax, double tmin,
			    double tmax,
			    double array1[][NYMAX])
			     //double** array1,double** array2,double** array3)

{

  int nxout = (int)(0.5+(xmax-xmin)/dxout)  +1;
  int ntout = (int)(0.5+(tmax-tmin)/dtout)  +1;
  if((nxout>NXMAX)||(ntout>NYMAX)){
    cerr <<"InOut::write_array2d:" 
         <<" nxout="<<nxout<<">=NXMAX="<<NXMAX<<" and/or ntout="<<ntout
	 <<" >=NYMAX="<<NYMAX<<" => Error!" << endl;
    exit(-1);
  }

  cout << "in write_array2d: "<<endl
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<", xmax="<<xmax <<endl
       <<" tmin="<<tmin<<", tmax="<<tmax <<endl
       <<" dxout="<<dxout<<", dtout="<<dtout <<endl
       <<" nxout="<<nxout<<", ntout="<<ntout
       <<endl;
  FILE  *outfile;                  
  //  filecheck (outfile,fname);
  outfile = fopen(fname,"w");

  fprintf(outfile, 
    "# x(km) \t t(min) \t rho(1/km) \t v(km/h) \t Q(1/h) \t a(m^2/s^2) \n");

  for (int it=0; it<ntout; it++)
  {
  
    for (int j=0; j<nxout; j++)
    {
      double x_km          = (xmin+(xmax-xmin)*j/nxout)/1000.;
      double t_h           = (tmin+(tmax-tmin)*it/ntout)/3600.;
      double rho_invkm     = array1[j][it];
      //double v_kmh         = array2[j][it];
      //double Q_invh        = array3[j][it];
      double v_kmh = 0;
      double Q_invh = 0;
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\t  %.6f\t  %f\n",
                      x_km, t_h, rho_invkm, v_kmh, Q_invh);
    }
    fprintf(outfile, "\n");        // newline mean new t step for gnuplot
  }

  fclose(outfile);

}




//#########################################################
//  special-purpose function  write_array2d (general purpose: above)
//#########################################################

  // definition of arrays for calling write_array2d: see HINT in .h file

void  InOut::write_array2d (char* fname, double dxout, double dtout, 
    double xmin, double xmax, double tmin, double tmax,
    double array1[][NYMAX], double array2[][NYMAX],
    double array3[][NYMAX], double array4[][NYMAX])
		      //double** array1,double** array2,double** array3)

{

  int nxout = (int)((xmax-xmin)/dxout);
  int ntout = (int)((tmax-tmin)/dtout);
  if((nxout>NXMAX)||(ntout>NYMAX)){
    cerr <<"InOut::write_array2d:" 
         <<" nxout="<<nxout<<">=NXMAX="<<NXMAX<<" and/or ntout="<<ntout
	 <<" >=NYMAX="<<NYMAX<<" => Error!" << endl;
    exit(-1);
  }

  cout << "in write_array2d: "<<endl
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<", xmax="<<xmax <<endl
       <<" tmin="<<tmin<<", tmax="<<tmax <<endl
       <<" dxout="<<dxout<<", dtout="<<dtout <<endl
       <<" nxout="<<nxout<<", ntout="<<ntout
       <<endl;
  FILE  *outfile;                  
  //  filecheck (outfile,fname);
  outfile = fopen(fname,"w");

  fprintf(outfile, 
    "# x(km) \t t(min) \t rho(1/km) \t v(km/h) \t Q(1/h) \t wMatrix\n");

  for (int it=0; it<=ntout; it++)
  {
  
    for (int j=0; j<=nxout; j++)
    {
      double x_km          = (xmin+(xmax-xmin)*j/nxout)/1000.;
      double t_h         = (tmin+(tmax-tmin)*it/ntout)/3600.;
      double rho_invkm     = array1[j][it];
      double v_kmh         = array2[j][it];
      double Q_invh        = array3[j][it];
      fprintf(outfile,  "%.6f\t  %.6f\t %.6f\t  %.6f\t  %.6f\t  %f\n",
                      x_km, t_h, rho_invkm, v_kmh, Q_invh, array4[j][it]);
    }
    fprintf(outfile, "\n");        // newline mean new t step for gnuplot
  }

  fclose(outfile);

} // InOut::write_array2d





//############################################################
//  function getvar
//############################################################

// interprets lines beginning with COMMENTCHAR as comments;
// reads the first word of other lines;
// inteprets the remaining contents also of these lines as comments
//!!! geht nicht mit referenz-Werten (F... File lesen)

void InOut::getvar(FILE *fp, double *pdouble)
{
  char   str[LINEMAX];               // Dummy string for comments 
  //  const char COMMENTCHAR = '%';
  if (fgets(str,LINEMAX,fp))
  {
    while ( ( (str[0]==COMMENTCHAR)||(str[0]==COMMENTCHAR2)||(str[0]=='\\'))
        && fgets(str,LINEMAX,fp)){
    }
    sscanf(str,"%lf", pdouble);
    printf("%.3f\t %s",*pdouble,str);
  }
} // InOut::getvar


void InOut::getvar(FILE *fp, int *pint)
{
  char   str[LINEMAX];    
  const char COMMENTCHAR = '%';

  if (fgets(str,LINEMAX,fp))
  {
    while ( ( (str[0]==COMMENTCHAR)||(str[0]==COMMENTCHAR2)||(str[0]=='\\'))
        && fgets(str,LINEMAX,fp)){
    }
    sscanf(str,"%i",pint);
    printf("%i\t %s",*pint,str);
  }
}


void InOut::getvar(FILE *fp, long *plong)
{
  char   str[LINEMAX];    
  const char COMMENTCHAR = '%';

  if (fgets(str,LINEMAX,fp))
  {
    while ( ( (str[0]==COMMENTCHAR)||(str[0]==COMMENTCHAR2)||(str[0]=='\\'))
        && fgets(str,LINEMAX,fp)){
    }
    sscanf(str,"%li",plong);
    printf("%li\t %s",*plong,str);
  }
}

void InOut::prepare_line(char line[]){
  // transform "csv" into spaces
  for (int j=0; j<LINEMAX; j++){
    if(line[j]==';'){
      line[j]=' ';
    }
  }
}

int InOut::is_not_white(char line[], int nchar)
{
  int not_white=false;
  for (int i=0; i<=nchar-1; i++)  // last always '\n'
  {
     // cout << " line["<<i<<"]="<<line[i] << "nonwhite = "
      //    <<( ( (line[i]!=' ') && (line[i] !='\n') )? 1 : 0) << endl;
    if( (line[i]!=' ') && (line[i] !='\n') && (line[i] !='\t') ) not_white=true;
  }
  return(not_white);
}

bool InOut:: is_data_line(char line[]){
  return  (line[0] != COMMENTCHAR) && 
    (line[0] != COMMENTCHAR2) && 
    (strlen(line)>0) ;
}




void  InOut::write_array2d (char* fname, double dxout, double dtout, 
    double xmin, double xmax, double tmin, double tmax,
    double array1[][NYMAX+1], double array2[][NYMAX+1],
			    double array3[][NYMAX+1], double array4[][NYMAX+1],
			    bool data3d)
			   

{

  int nxout = (int)((xmax-xmin)/dxout)+1; // write_array2d (multivariate data)
  int ntout = (int)((tmax-tmin)/dtout)+1; // write_array2d (multivariate data)
  if((nxout>NXMAX)||(ntout>NYMAX)){
    cerr <<"InOut::write_array2d:" 
         <<" nxout="<<nxout<<">=NXMAX="<<NXMAX<<" and/or ntout="<<ntout
	 <<" >=NYMAX="<<NYMAX<<" => Error!" << endl;
    exit(-1);
  }

  cout << "in write_array2d: "<<endl
       << "fname = "<<fname <<endl
       <<" xmin="<<xmin<<", xmax="<<xmax <<endl
       <<" tmin="<<tmin<<", tmax="<<tmax <<endl
       <<" dxout="<<dxout<<", dtout="<<dtout <<endl
       <<" nxout="<<nxout<<", ntout="<<ntout
       <<endl;
  FILE  *outfile;                  
  //  filecheck (outfile,fname);
  outfile = fopen(fname,"w");

  if(data3d){
    fprintf(outfile, "# x(inp units) \t t(inp units) \t rho \t v(inp units) \t Q\t wMatrix\n");
  }
  else{
      fprintf(outfile, "# x(km) \t t(h) \t\t rho(1/km) \t v(km/h) \t Q(1/h) \t wMatrix\n");
  }
  
  // (martin nov07) no longer case distinction necessary
  
  for (int it=0; it<ntout; it++){
     for (int j=0; j<nxout; j++) {
       double x          = (xmin+(xmax-xmin)*j/(nxout-1));
       double t        = (tmin+(tmax-tmin)*it/(ntout-1));
       double rho    = array1[j][it];
       double v         = array2[j][it];
       double Q       = array3[j][it];
       fprintf(outfile,  "%.6f\t  %.6f\t %.6f\t  %.6f\t  %.6f\t  %f\n",
                      x, t, rho, v, Q, array4[j][it]);
     }

     fprintf(outfile, "\n");        // newline mean new t step for gnuplot
  }

  fclose(outfile);

}


