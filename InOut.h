#ifndef IN_OUT_H
#define IN_OUT_H

//(jun 04) ACHTUNG: InOut gibt nun n_array als Zahl der Eintraege raus!
// max. Index also n_array-1  !!


/// Allgemeine I/O-Routinen
class InOut{

 public:

  InOut();
  static const int NXMAX=1000;
  static const int NYMAX=101;
  static const int NTMAX=1441;

  static const int LINEMAX=1500;
  //static const int MAXSTR=256;
  static const char COMMENTCHAR = '%';
  static const char COMMENTCHAR2 = '#';

  int getNumberOfLines(char* fname);
  void  get_array (char* fname, int & n_array, double array1[]);
  void  get_array (char* fname, int & n_array, int array1[]);
  void  get_array (char* fname, int & n_array, 
  		   double array1[], double array2[]);
  void  get_array (char* fname, int & n_array, 
		   double array1[], double array2[], double array3[]);

  void  get_array (char* fname, int & n_array, 
		   int array1[], double array2[], double array3[]);

  void  get_col (const char* fname, int col,  int & n_array, int array1[]);
  void  get_col (const char* fname, int col,  int & n_array, double array1[]);
  void  getChar_col (const char* fname, int col,  int & n_array, char array1[]);

  void  write_array (char* fname, int nData, const double data_col1[],
                      char* name_col1);
  void  write_array (char* fname, int nData, const double data_col1[],
                     const double data_col2[], 
                     char* name_col1, char* name_col2);
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],
                     char* name_col1, char* name_col2, char* name_col3);
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],
                     char* name_col1, char* name_col2, char* name_col3,
		     char* commentline);
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4);
  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8);

  void  write_array (char* fname, int nData, 
                     const double data_col1[],const double data_col2[], 
                     const double data_col3[],const double data_col4[],
                     const double data_col5[],const double data_col6[],
                     const double data_col7[],const double data_col8[],
                     char* name_col1, char* name_col2, 
                     char* name_col3, char* name_col4, 
                     char* name_col5, char* name_col6, 
                     char* name_col7, char* name_col8,
		     char* commentline);


  void  write_array2d (char* fname, int nt, int nx, double dt, double dx,
  double array1[NXMAX+1][NTMAX+1],
  double array2[NXMAX+1][NTMAX+1],
  double array3[NXMAX+1][NTMAX+1]);

  void  write_array2d (char* fname, double dxout, double dtout, 
                       double xmin, double xmax, 
                       double tmin, double tmax, 
                       double array1[NXMAX+1][NTMAX+1]);


  void  write_array2d (char* fname, double dxout, double dtout, 
    double xmin, double xmax, double tmin, double tmax,
                       double array1[][NTMAX+1],
                       double array2[][NTMAX+1],
                       double array3[][NTMAX+1],
		       double array4[][NTMAX+1]);


  void  get_array2d_col (const char* fname, int col, 
		       double& xmin, double& dx, int& nx, 
		       double& ymin, double& dy, int& ny, 
		       double array[NXMAX][NYMAX]);


  static void getvar(FILE *fp, int* pint);
  static void getvar(FILE *fp, long* plong);
  static void getvar(FILE *fp, double* pdouble);

 private:

  int is_not_white(char line[], int nchar);
  void prepare_line(char line[]);
  bool is_data_line(char line[]);


};

#endif // IN_OUT_H
