#ifndef FPE_H
#define FPE_H
class FPE{

 public:

  FPE();
  FPE(double dx, double xWidth,int choice_method);

  void get_modelparams(const char fname[]);
  void calc_rhs(const double f[],
              double F1[],
              double S1[],
		bool downwind_diff);

 private:

  char fName[256];

  double A;
  double B;

  double xmin;
  double xWidth;
  double dx;
  int nx;

};

#endif // FPE_H
