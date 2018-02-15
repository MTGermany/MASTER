#ifndef KKL_H
#define KKL_H
class KKLparams
{
public:
  KKLparams();
  void getKKLparams (char namepar[]);
  double c0;
  double mu;
  double tau;  
  double lee_E; 
  int lee_theta;
  double kerner_rhoi;
  double kerner_b;
};

#endif // KKL_H
