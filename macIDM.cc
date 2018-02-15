
// other parameters from .inp file: v0, rhomax, t, a=v0/tau

IDMparams::IDMparams()
{
  b  = 2.;         // everyday braking deceleration
  s0 = 0.;          // minimal gap (used also in calc_rhs of GKT model)
  delta = 4;       // exponent
}


void IDMparams::getIDMparams(char namepar[])
{
  FILE * fp;                    
  char   param_file[MAXSTR];   
  sprintf(param_file,"%s.IDM",namepar);
  cout << "\nIDM model chosen:" << endl
       << "Reading params of IDM model from file "<< param_file
       << endl;
  fp = fopen(param_file,"r");
  filecheck(fp, param_file);

  getvar(fp,&b);
  getvar(fp,&s0);
  getvar(fp,&delta);

  fclose(fp);
}

void   IDMparams::calc_eq( double vwtab[],  double veqtab[])

    // Calculates equilibrium velocity of GKT and IDM with finite s0
    // and free-acc exponent delta
    // uses numeric iteration procedure and
    // !! calculates THE WHOLE FIELD veq
    // NO resignation effect respected up to now
{

  int    ir;
  int    itmax       = 100;
  double dtmax       = 2;       // dtmax,dtmin (in s) v dependent time steps
  double dtmin       = 0.01;    // for relaxation method

  veqtab[0]     = v0;
  double v_it;

  for(ir=1; ir<=NRHO; ir++)
  {
      double rho    = rhomax*ir/NRHO;
      double s      = 1./rho - 1./rhomax;
      double v0rho  = vwtab[ir];  
      v_it          = veqtab[ir-1];

      for (int it=1; it<=itmax; it++){
        double dtloc    = dtmax*v_it/v0rho + dtmin;
        double dxstar = (s0 + Tr * v_it)
                        + s0*sqrt((v_it+0.000001)/v0rho);
	// !! Acceleration  propto v0rho/tau0, braking propto v0/tau0
        double a      = (s>s0) 
            ? v0rho/tau0 * (1-pow(v_it/v0rho, delta)) - v0/tau0*SQR(dxstar/s)
            :0.;
        v_it         += dtloc*a;
        if((v_it<0)||(s<s0)) v_it=0;

        if(false) 
          cout << "rho="<< rho<< " s="<<s<<" dxstar="<<dxstar
               <<" v0rho="<<v0rho<<" a="<<a<<" v_it="<<v_it<<endl;
	}
      veqtab[ir] = v_it;
  }
   // cout << "denom="<<denom<<" vw="<<vw
  //        <<" vtilde="<<vtilde<<" veq="<<veq<<endl; 
}
