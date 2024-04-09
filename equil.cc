
double tiny_value=1.e-6;

//#######################################################
// minsmooth and maxsmooth (needed for general FD of LWR model, choice_model=9)
//#######################################################

double minsmooth(double x1, double x2, double dx){
  return 0.5*(x1+x2) - sqrt(0.25*(x1-x2)*(x1-x2) + dx*dx);
}

double maxsmooth(double x1, double x2, double dx){
  return 0.5*(x1+x2) + sqrt(0.25*(x1-x2)*(x1-x2) + dx*dx);
}




// #######################################################
void calc_tables(char namepar[])
// #######################################################

  /* 
  (1) Calculates some tabulated functions
     + rho dependend functions:
       - veqtab   = equilibrium velocity for homogeneous traffic
       - sqrtAtab = sqrt of the A(rho) function
       - v0factab = v0(rho)/v0(0) in the KKL model
       - high_denstab (only GKT with resignation factor or 
         high density corection 
     + Maximum equilibrium flow Qmax and corresponding abszissa rhoQmax
     + flow dependend functions:
       rho_freetab = free branch and rho_congtab = congested branch
       of homogeneous density for given flow
     + Btab = velocity-difference dependent GGKT interaction factor 
       B(delta_v)

  (2) Writes rho, A(rho), veq(rho), and Qeq(rho) in file <namepar>.tab


   choice_model=0:  GKT
   choice_model=1:  reduced desired velocity and additional diffusion
                    for high densities; NO decrease of interaction).
   choice_model=2:  Kerner-Konhaeuser model
   choice_model=3:  macroscopic IDM
   choice_model=4:  macroscopic OVM (with ve as Bando and max(v0,s/T))
   choice_model=5:  Black-Scholes eq as example for 1-field parab. eq
   choice_model=6:  GKT with resignation factor
   choice_model=7:  GKT with jam distance s0 and free-acc exponent delta=4
   choice_model=8:  Macro-3-Phase model (mar08)
   choice_model=9:  LWR with parabolic FD (apr16)
   choice_model=10: SGM (2018) (speed gradient model Jian 2002) (=> ISTTT23)
   choice_model=21: VMM model (2017)
   choice_model=50: FPE
   choice_model=51: FPE innov (see calc_rhs.cc)



  **********************************************************************/


{

  FILE * fptab;                  
  char   tab_file[MAXSTR];       
  int ir;

  // ########################################################
  // First rho loop: Calculate v0factab (KKL model, SGM) and high_denstab 
  // (GKT variants, not used for original GKT model)
  // v0factab in [0,1] lowers desired velocity;
  // high_denstab used EXCLUSIVELY for diffusion; 
  // ########################################################

  if(true) cout << "\nin calc_tables of equil: First loop" << endl;

  if(choice_model==10){// SGM independent!
    sgm.calc_eq(veqtab); // => veqtab[ir] (!! veqtab here and in SGM class
  }

  else{ // all other models connected=>change!

   for(ir=0; ir<=NRHO; ir++)
   {
    double rho = rhomax*ir/NRHO;

    if (choice_model==1)   // GKT with high-density correction          
    {
      double rho_switch      = 0.75*rhomax;        //original 0.75
      double drho            = 0.05*rhomax;        //original 0.05
      high_denstab[ir]       = 1./(1.+ exp(-(rho-rho_switch)/drho ));
      v0factab[ir]           = 1.-high_denstab[ir];
    }

    else if (choice_model==6)   // GKT with "resignation factor"  
    {                                          
                                               // A5               A9
      double rhoc_v0      = 0.48*rhomax;        // 0.4              0.42
      double drho_v0      = 0.12*rhomax;       // 0.08             0.12
      double fac_resign   = 0.1;               // 0.2   (in units of v0)

      
      // addtl. high dens corr bringt Verschlechterung -> nix!
      
      v0factab[ir]      = 1.-(1-fac_resign)/(1.+ exp(-(rho-rhoc_v0)/drho_v0 ));
      double rhoc_diff    = 0.4*rhomax;   // 0.4; this and
                              //  following only needed if  diffusion is on
      double drho_diff    = 0.08*rhomax;      
      high_denstab[ir]  = 1./(1.+ exp(-(rho-rhoc_diff)/drho_diff ));
    }

    else if (choice_model==2)  // KKL model 
    {
      high_denstab[ir]= 0.;   // no additional diffusion to that of KK-model

      // Comment out either the 4 lines below (1) -> Ve(rho) of Lee
      // or the 2 lines below (2) -> Ve(rho) of Kerner-Konhaeuser

      /*
      // (1) original KK equilibrium relation

      double d         = 1./(1.+exp((1. -kkl.kerner_rhoi/rhomax )
                           /kkl.kerner_b));
      v0factab[ir]     = 1./(1.+exp((rho/rhomax-kkl.kerner_rhoi/rhomax)
                           /kkl.kerner_b)) - d;
      */

      // (2) or alternatively equilibrium relation of 
      // Lee/Lee/Kim, PRL81,1130 (1998) 
      // Comment out one of these!

      v0factab[ir]     = (1. - rho/rhomax) 
             / (1 + kkl.lee_E * pow( rho/rhomax, kkl.lee_theta));

    }


    else   // Other models including original GKT and macroscopic IDM and OVM
    {
      high_denstab[ir]=0.;
      v0factab[ir] = 1.;     

    }
   }
  }

  // ********************************************************************
  // Second loop: Calculate variance prefactor A(rho): theta=A(rho) V^2
  // ********************************************************************

  cout << "in calc_tables: Second  loop GKT-A(rho)" << endl;
  if(choice_model==10){cout<<"SGM: do nothing"<<endl;}



     // const A

  else if      ((choice_A == 0) || (choice_A==1)) { 
    for(ir=0; ir<=NRHO; ir++)  sqrtAtab  [ir] = sqrt(A0);
    Arhomax = A0;      //A(rhomax); needed in velocity eq of GKT formula
  }

  else if (choice_A == 2)   
  {
    char   infname[MAXSTR];   // name of the file with the A parameters
    double pos_rel;
    double width_rel;      // the larger, the less sharp is the step
    FILE  *infp;
    sprintf(infname,"%s.A",namepar);

    if( !( infp=fopen(infname,"r")) )  
    {
      cout << "No file "<< infname << " detected -> "
           << "Using default parameters"<< endl << "for A function:"
           << "pos_rel = 0.27;  width_rel = 0.1 " << endl;
      pos_rel   = 0.27;
      width_rel = 0.1;
    }

    else
    {
      cout << endl << "Reading parameters of A function from file " 
           << infname << endl;
      getvar(infp,&pos_rel);  
      getvar(infp,&width_rel);
      fclose(infp);
    }

    for(ir=0; ir<=NRHO; ir++)
    {
      double rho      = rhomax*ir/NRHO;
      double arg      = (rho/rhomax - pos_rel)/width_rel; 
      sqrtAtab  [ir]  = sqrt(A0 + dA * ( tanh(arg)+1.)); 
    }

    Arhomax = A0 + 2.*dA; 
  }



  // **********************************************************
  // Third rho loop: Calculation of equil.  values (with v0 from .inp!)
  // **********************************************************

  if(choice_model==10){cout<<"SGM: Do nothing in third loop"<<endl;}


  else if ((choice_model==0) ||   // original GKT
      (choice_model==1) ||   // GKT with high-density corr
      (choice_model==6))     // GKT with resignatione effects
  {
    for( ir=0; ir<=NRHO; ir++)
    {
      double rho           = rhomax*ir/NRHO;
      calc_eq_GKT(rho,  v0*v0factab[ir], SQR(sqrtAtab[ir]), veqtab[ir]);
    }
  }  

  else if(choice_model==9){// LWR (with transition triang-parabolic FD controlled by A0)
    for (ir=0; ir<=NRHO; ir++){
      double rho           = rhomax*ir/NRHO;
      double veParabola=v0*(1.-rho/rhomax);
      double veTriang
	=max(0.,minsmooth(v0, (1/max(rho,tiny_value)-1/rhomax)/Tr, 0.8*veParabola*A0));
      veqtab[ir]=A0*veParabola+(1-A0)*veTriang;
    }
  }
 
  else if (choice_model==2)  // KKK model 
    for( ir=0; ir<=NRHO; ir++)
      veqtab[ir]=v0 * v0factab[ir];


  else if ((choice_model==3)   // "macroscopic" IDM 
	|| (choice_model==7))  // GKT with s0 and accel exponent delta
  {
    double vwtab [NRHO+1];
    for( ir=0; ir<=NRHO; ir++) vwtab[ir] = v0*v0factab[ir];
    idm.calc_eq(vwtab, veqtab);
  }


  else if (choice_model==4) // "macroscopic" OVM
  {
    bool useMaxV0_sdurchT=true; //!!!
    //bool useMaxV0_sdurchT=false;

    if(!useMaxV0_sdurchT){ // use original ve(rho) relation of Bando et. al.
      veqtab[0] = ovm.V0 * ( tanh(ovm.c1)+1 );
      for( ir=1; ir<=NRHO; ir++){
        double rho=rhomax*ir/NRHO;
        veqtab[ir] = ovm.V0
	  * (tanh(ovm.c1) + tanh(ovm.c0*(1/rho-1/rhomax)-ovm.c1))
	  / (tanh(ovm.c1)+1);
        if (veqtab[ir] <= 0.0) veqtab[ir] = 0.0;
      }
    }
    else{               // use ve(rho)=max(v0,s/T)
      //double v0=35;  //!!! useMaxV0_sdurchT: v0 from .inp
      double T =1.5; //!!!  useMaxV0_sdurchT: OVM T ad-hoc!
      // double s0=0; //!!! OVM s0 ad-hoc!
      double rel_capacDrop=0.2;  // Qdrop/Qmax //0  0.2
      double rhoc=(rel_capacDrop+1)/(T*v0+1./rhomax);
      veqtab[0] = v0;
      for( ir=1; ir<=NRHO; ir++){
        double rho=rhomax*ir/NRHO;
	veqtab[ir] = (rho<=rhoc) ?  v0 : (1/rho-1/rhomax)/T;
      }
    }

  }  
  else if ((choice_model==5)||(choice_model==21)||(choice_model==50)||(choice_model==51)){ 
  // Black Scholes, FPE and other models: nothing needed?
  }

  else if (choice_model==8){
    veqtab[0] = mac3phases.get_v0();
     for(int ir=1; ir<=NRHO; ir++){
       cout <<" ir="<<ir<<endl;
        double rho=rhomax*ir/NRHO;
        veqtab[ir] = mac3phases.ve_rho(rho);
      }
  }

  
  else error(" equil.cc: sorry, choice_model<0 or neq 21,50,51 or >10 not implemented");

  if (test_timestep) {
    for( ir=0; ir<=NRHO; ir++)
      if ((ir<=5) || (NRHO-ir <=5))
        cout<<"ir = "<<ir
            <<" ve = "<<veqtab[ir] 
            << endl;
  }

  // ********************************************
  // Calculate maximum of fundamental diagram:
  // Qmax and abszissa rhoQmax
  // ********************************************

  Qmax=-1.;
  ir=1;

  while(veqtab[ir] * rhomax*ir/NRHO > Qmax)
  {
    Qmax = veqtab[ir] * rhomax*ir/NRHO;
    ir++;
  }
  rhoQmax = rhomax*ir/NRHO;
  cout << "calc_tables: Qmax="<<Qmax<<" rhoQmax="<<rhoQmax<<endl;

  /* *********************************************************
   Forth rho loop: Calculation of  the inverse of the Qe(rho)
   and ve(rho) relations:
   tables rho_freetab[iQ] = rho^-1(Qe) if Q<=Qmax
          rho_freetab[iQ] = rhoQmax    if Q>Qmax
   equidistant Q values from 0..Qtabmax>Qmax
   rho_congtab analog for congested branch
   *********************************************** */

  double Qtabmax   = QFACTAB*Qmax;     // max. Q value of tables (>Qmax!)
  int    irhoQmax   = (int) (NRHO*rhoQmax/rhomax);
  double ord_v_tab       [NRHO+1];  
  double ord_rho_freetab [NRHO+1];  
  double ord_rho_congtab [NRHO+1];  
  double veqtab_back     [NRHO+1];  
  double Qfreetab        [NRHO+1];  
  double Qcongtab        [NRHO+10];  

  // make tables for the intpextp routine
  // The abszissa table must be monotonuously increasing !!!

  for (int ir1=0; ir1<=NRHO; ir1++) {
    int ir = NRHO-ir1;
    ord_v_tab   [ir1] = rhomax*ir/NRHO;         // ordinate rho(v)
    veqtab_back [ir1] = veqtab[ir];             // increasing abszissa v
    //cout
    //    <<"veqtab_back [ir1]=" << veqtab_back [ir1] 
    //    <<" ord_v_tab [ir1]="<<ord_v_tab [ir1]
    //    << endl;
  }

  for (int ir=0; ir<=irhoQmax; ir++) {
    ord_rho_freetab[ir] = rhomax*ir/NRHO;
    Qfreetab   [ir] = ord_rho_freetab[ir] * veqtab[ir];
  }

  for (int ir1=0; ir1<=NRHO-irhoQmax; ir1++) {
    int ir = NRHO-ir1;
    ord_rho_congtab [ir1] = rhomax*ir/NRHO; 
    Qcongtab        [ir1] = ord_rho_congtab[ir1] * veqtab[ir];
  }

  // intpextp makes automatically constant part for Q>Qmax

  for (int iv=0; iv<=NRHO; iv++) 
    rho_vtab [iv] = 
      intpextp(veqtab_back, ord_v_tab,  NRHO+1, v0*iv/NRHO);

  for (int iQ=0; iQ<=NRHO; iQ++) 
    rho_freetab[iQ] = 
      intpextp(Qfreetab, ord_rho_freetab,  irhoQmax+1, Qtabmax*iQ/NRHO);


  for (int iQ=0; iQ<=NRHO; iQ++) 
    rho_congtab[iQ] = 
      intpextp(Qcongtab, ord_rho_congtab, NRHO-irhoQmax+1, Qtabmax*iQ/NRHO);

  //if(test_timestep){
  if(true) {

    cout << "calc_tables after 4th loop: "<<endl<<endl;

    for (int iv=0; iv<=10; iv++) {
    cout << " v="<<v0*iv/NRHO
         << " rho_vtab[iv]="<<rho_vtab[iv]<<endl;
    }


    for (int iQ=0; iQ<=NRHO; iQ++) {
      cout << " Q="<<Qtabmax*iQ/NRHO
         << " rho_freetab[iQ]="<<rho_freetab[iQ]
         << " rho_congtab[iQ]="<<rho_congtab[iQ]
         << endl;
    }
    
  }


  // ***********************************************************
  // last rho loop: Output in file
  // ************************************************************

  sprintf(tab_file,"%s.tab",namepar);
  fptab=fopen(tab_file,"w");

  // write title line ( "#" = comment for gnuplot)

  Q_max=0.0;

  fprintf(fptab, "# rho\t    A\t        veq\t       Qeq\t  \n");

  for(ir=0; ir<=NRHO; ir++)
  {
    double rho              = rhomax*ir/NRHO;
    if (rho*veqtab[ir] > Q_max) Q_max = rho*veqtab[ir];
    fprintf(fptab, "%.1f\t  %.4f\t %.2f\t       %.1f \n",
       1000.*rho, SQR( sqrtAtab[ir]), 
       3.6*veqtab[ir], 3600.*rho*veqtab[ir]);
  }

  if (test_timestep)
    cout << "rho, SQR(sqrtAtab), veqtab, Qeqtab "
         << "saved in " << tab_file << endl;
  fclose(fptab);


  /**********************************************************************
   (2.) Calculation of the velocity-difference dependent interaction factor
        B(delta_V); 
        B(delta_V) defined such that B(0) = 1
  **********************************************************************/

  double delta_dv = (DVMAX-DVMIN)/NDV;
    
  for(int i=0; i<=NDV; i++)
    {
    double dv       = DVMIN + i*delta_dv;
    double erf      = ( (dv<XERFMAX) && (dv>XERFMIN) )
                    ? intp(erftab, NERFTAB+1, dv, XERFMIN, XERFMAX)
                    : (dv<0) ? 0 : 1;
    double gauss    = 1./sqrt(2.*PI)*exp(-SQR(dv)/2.);
    Btab[i] = 2.* ( dv*gauss + (1.+SQR(dv))*erf);
    //    if(test_timestep) cout << "dv = " << dv << "erf = " << erf 
    //   << "Btab= " << Btab[i] << endl;
    }

  }


/********************************************************************/
/*  function  calc_eq_GKT                                              */
/********************************************************************/

void calc_eq_GKT(double rho, double vw, double A,  double &veq)

    // Calculates equilibrium velocity of GKT and IDM (the same for s0=0)
    // directly with the solution of 
    // the quadratic equation
    // In case of space-dependent model aprameters use that of .inp

{
  if(rho<1/HUGE_VALUE) rho=1/HUGE_VALUE;
  if(rho>rhomax-1/HUGE_VALUE) rho=rhomax-1/HUGE_VALUE;

  double denom   =  1 - rho/rhomax;
  double Arhomax_loc = SQR( (sqrtAtab[NRHO]));
  double vtilde  =  denom / (rho*Tr) * sqrt(Arhomax_loc/A); 
  veq            = SQR(vtilde)/( 2. * vw) 
                     * (-1. + sqrt(1. + 4.*SQR(vw/vtilde)));
  // cout << "denom="<<denom<<" vw="<<vw<<" vtilde="<<vtilde<<" veq="<<veq<<endl; 
}




