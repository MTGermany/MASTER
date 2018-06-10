




/********************************************************************/
/*  function calc_rhs                                               */
/********************************************************************/


void calc_rhs (int choice_model, bool downwind_diff,
              double rho[],double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[], 
              bool show_calc_rhs)

  /* Calculates the right-hand sides (fluxes Fi and sources Si)
     of macroscopic traffic models of the form
              ________________________________________________
              |  d_t rho(x,t)  = -d_x F1 + S1 (ramp flows)    |
              |  d_t Q(x,t)    = -d_x F2 + S2                 |
              ------------------------------------------------

     Input: 
       - choice_model: defines the traffic model; choice_model==0: GKT model,
         see equil.cc for details!
       - downwind_diff: whether downwind differences should be used
         for calculating derivatives (only relevant if fluxes with gradients)
       - The fields rho,Q. They will be modified only, if 
         something unphysical happens
       - switch "show_calc_rhs" to show some debugging information

     Output: 
       - The fluxes F1,F2, and the sources S1,S2 of the model.
       - The sources include ramps: S1 = ramp flow/merging length 
                                    S2+= S1*(velocity when merging) 
       - Flow-conserving bottlenecks realized by gradients of the
         parameters T (=Tr) or v0, given by the (global) arrays
         Tr_loc and v0_loc
   */


{

  int    i;
  int    nx   =(int)(xmax/dx);
  double v[NXMAX+1];                         // velocity
  double rhodelta[NXMAX+1];                  // nonlocal density
  double vdelta[NXMAX+1];                    // nonlocal velocity


  // ######################################################################
  // Standard GKT model
  // d_t rho = - d_x Q + ramp source terms,
  // d_t Q   = - d_x (rho A V^2) + (rho V_0 - Q)/tau0
  //           - v0 rho' T^2 A(rho) Q^2 B(dvrel) 
  //             / [tau0 A(rhomax) (1-rho'/rhomax)^2]      
  //
  // with A=A(rho) variance prefactor, 
  // rho'     = rho at anticipated position, (V' analog)
  // B(dvrel) = gas-kinetic "Boltzmann" interaction prefactor [B(0)=1]
  // dvrel    = velocity difference (V - V') in units of the sqrt of the
  //            velocity variance A V^2

  // ######################################################################

  bool testIsolatedGKT=false;

  if  ( testIsolatedGKT&&(!trucks.varTruck) && (choice_model==0) ){
    calc_rhsGKT (choice_model, downwind_diff, rho,Q,F1,F2,S1,S2, show_calc_rhs);
  }


  else if  ( (!trucks.varTruck) && (choice_model==0) ) {

    // Check input (rho>0, Q>=0) and calculate V

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<TINY_VALUE)        rho[i] = TINY_VALUE;
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]   =rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }

    // calculate advanced rho and v fields
    
    shift_in_x (antic_factor, v, rho, rhodelta);
    shift_in_x (antic_factor, v, v,   vdelta);

    // calculate rhs of flow equation

    for (i=0; i<=nx; i++)
    {
      double sqrtA     = intp(sqrtAtab, NRHO+1, rho[i], 0, rhomax);
      double denom     = 1-rhodelta[i]/rhomax;
      double dvrel     = (v[i]-vdelta[i]) / (sqrtA * sqrt2 * v[i]  );   
      double bolzm_fact= ( (dvrel>DVMIN) && (dvrel<DVMAX) )
                       ? intp( Btab, NDV+1, dvrel, DVMIN, DVMAX)
                       : (dvrel>0) ? 2.*SQR(dvrel) : 0.; 

      // calculate F1, F2 and S2; S1 calculated in ramp-source section
      // below

      F1[i] = Q[i];
      F2[i] = rho[i] * SQR(v[i]) * (1.+ SQR(sqrtA)); // hat wenig Einfluss
      //F2[i] = rho[i] * SQR(v[i]) * (1.);   //martin08
      double S2free = (rho[i]*v0_loc[i] - Q[i]) / tau0;
      //double S2brake= v0_loc[i] * rhodelta[i] //orig
      double S2brake= v0_loc[i] * rho[i]  //martin08
	//* SQR( Tr_loc[i]*sqrtA*Q[i]) * bolzm_fact //orig
                      * SQR( Tr_loc[i]*sqrtA*rhodelta[i]*v[i]) * bolzm_fact //martin08
                      / (SQR(denom)*tau0*Arhomax);
      S2[i] = S2free  - S2brake;
    }
  }


  // #######################################################
  // GKT model with high-density correctionn 
  // #######################################################

  else if  ( (!trucks.varTruck) && (choice_model==1) )
  {

    // max. diffusion (m^2/s) for rho=rhomax
    // set to 0 if more stable flow-violating
    // method in timestep is used

    double Dplusmax=200;  //!!!
    const double durch_dx = 1./dx;

    // Check input (rho>0, Q>=0) and calculate V

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<TINY_VALUE)        rho[i] = TINY_VALUE;
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]   =rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }

    // calculate advanced rho and v fields
    
    shift_in_x (antic_factor, v, rho, rhodelta);
    shift_in_x (antic_factor, v, v,   vdelta);

    // calculate rhs of flow equation

    for (i=0; i<=nx; i++)
    {
      double sqrtA     = intp(sqrtAtab, NRHO+1, rho[i], 0, rhomax);
      double denom     = 1-rhodelta[i]/rhomax;
      double dvrel     = (v[i]-vdelta[i]) / (sqrtA * sqrt2 * v[i]  );   
      double bolzm_fact= ( (dvrel>DVMIN) && (dvrel<DVMAX) )
                       ? intp( Btab, NDV+1, dvrel, DVMIN, DVMAX)
                       : (dvrel>0) ? 2.*SQR(dvrel) : 0.; 

      // calculate local high-density correction

      double corr_factor = intp(high_denstab, NRHO+1, rho[i],0,rhomax);
      double v0_corr     = v0_loc[i] * (1-corr_factor); // reduced v0
      double D_corr      = Dplusmax * corr_factor;      // additl diffusion
 
      // calculate F1, F2 and S2; S1 calculated in ramp-source section
      // below

      F1[i] = Q[i];
      F2[i] = rho[i] * SQR(v[i]) * (1.+ SQR(sqrtA));
      double S2free = (rho[i]*v0_corr - Q[i]) / tau0;
      double S2brake= v0_corr  * rhodelta[i]
                      * SQR( Tr_loc[i]*sqrtA*Q[i]) * bolzm_fact
                      / (SQR(denom)*tau0*Arhomax);
      S2[i] = S2free  - S2brake;

    // include high-density diffusion

      if (downwind_diff){
	 F1[i] += - D_corr * ( rho[i+1] - rho[i]) * durch_dx;
         F2[i] += - D_corr * ( Q[i+1]   - Q[i]  ) * durch_dx;
      }

      else{
	 F1[i] += - D_corr * ( rho[i] - rho[i-1]) * durch_dx;
         F2[i] += - D_corr * ( Q[i]   - Q[i-1]  ) * durch_dx;
      } 
    }
  }

  // ****** end GKT model with high-dens correctionn *****



  // **************************************************************
  // Other variants of the GKT model: With time-dependent truck percentage
  // (trucks.varTruck is tue), 
  // High-density corection (choice_model==1), 
  // Resignation factor (choice_model==6) (controlled in equil.cc), or
  // with s0 and free-acceleration exponent delta (choice_model==7).
  // **************************************************************

  else if  ( ((trucks.varTruck)&&(choice_model==0))
          || ((trucks.varTruck)&&(choice_model==1)) 
          || (choice_model==6) || (choice_model==7))
  {


    // Test input and calculate v

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<TINY_VALUE)  rho[i] = TINY_VALUE;
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]=rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }


    // calculate advanced rho and v fields
    
    shift_in_x (antic_factor,  v,rho,  rhodelta);
    shift_in_x (antic_factor,  v,v,    vdelta);
    //shift_integro (antic_factor,  v,rho,  rhodelta);
    //shift_integro (antic_factor,  v,v,    vdelta);


    for (i=0; i<=nx; i++)
    {
      double v0fac     = intp(v0factab, NRHO+1, rhodelta[i],0,rhomax); 
      double vw        = v0fac * v0_loc[i] * trucks.v0TruckCorr;
      double sqrtA     = intp(sqrtAtab, NRHO+1, rho[i],     0,rhomax);

      double denom     = 1-rhodelta[i]/(rhomax*trucks.rhomaxTruckCorr);
      double Csqrv     = vw*rhodelta[i]
                         * SQR((idm.s0+v[i]*Tr_loc[i]*trucks.TrTruckCorr)/denom) 
                         / (tau0*Arhomax);
      double dvrel     = (v[i]-vdelta[i]) 
                   / (sqrtA * sqrt2 * v[i]  );   // simplified dvrel 
      //  / (sqrtA *sqrt( SQR(v[i])+SQR(vdelta[i]) ) );  // original dvrel

      double bolzm_fact = ( (dvrel>DVMIN) && (dvrel<DVMAX) )
                       ? intp( Btab, NDV+1, dvrel, DVMIN, DVMAX)
                       : (dvrel>0) ? 2.*SQR(dvrel) : 0.; 

      // output of calc_rhs

      F1[i] = Q[i];
      F2[i] = rho[i] * SQR(v[i]) * (1.+ SQR(sqrtA));
      double S2free = (choice_model!=7)  
        ? (rho[i]*vw - Q[i]) / tau0
	: rho[i]*vw/ tau0 * (1.-SQR(SQR(Q[i]/(rho[i]*vw))));
      S2[i] = S2free  - Csqrv * SQR(rho[i]*sqrtA ) * bolzm_fact;   
    }

    // additional diffusion for high-density correction
    //!! switch off, if you want more stable but conservation-violating
    //!! method in "timestep.cc" (search for "Dplusmax")
    
    if((choice_model==1) || (choice_model==6))
    {
      // very sensitive; 100 does not work, 500 too much       
      // damping, damping only Q does not work

      //double Dplusmax=0;  
      double Dplusmax=200;  

      int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                      // ishift=1 for downwind diffrences

      for (i=1; i<=nx-1; i++) // no diffusion at the boundaries 
      {
         double Dactual = Dplusmax*(intp(high_denstab, NRHO+1, rho[i],0,rhomax));
	 F1[i] += - Dactual*( rho[i+ishift]-rho[i+ishift-1])/dx;
         F2[i] += - Dactual*( Q[i+ishift]-Q[i+ishift-1])/dx;
	 //  double vminVirt=-0.2;             // virtual desired veloc at v=0
	 // S2[i] += rho[i]*vminVirt/tau0; 
      }
    }
    
  }

  // ************* end Variants for GKT model *********************




  /***************************************************************
       Kerner-Konhaeuser model in flux representation:
 
        d_t Q = - d_x[ rho(v^2+c0^2) - mu d_x v] + (rho VKK-Q) / tau
        KK: tau=5 s, rhomax=175/km, vf=v0=100.8 km/h -> see "KK" in equil.c
        Lee: tau=30 s, rhomax = 140 veh/km, v0=(120/3.6) m/s
  ***************************************************************/

  else if (choice_model==2)  
  {
    for (i=0; i<=nx; i++){
      if (rho[i]<TINY_VALUE)  rho[i] = TINY_VALUE;
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]=rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }


    int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                        // ishift=1 for downwind diffrences

    double flux_diffus[NXMAX+1];  // diffusion comes into the flux F2
    flux_diffus[0]=flux_diffus[nx]=0; // no diffusion at the boundaries 
    for (i=1; i<=nx-1; i++){
       flux_diffus[i] = - kkl.mu*( v[i+ishift]-v[i+ishift-1])/dx;
    }

    for (i=0; i<=nx; i++)
    {
      // v0fac = ve(rho)/ve(0) -> v0factab calculated in equil.cc
      double v0fac     = intp(v0factab, NRHO+1, rho[i],0,rhomax); 
      double vw        = v0fac * v0_loc[i]* trucks.v0TruckCorr;
      F1[i] = Q[i];
      F2[i] = rho[i] * ( SQR(v[i]) + SQR(kkl.c0) ) + flux_diffus[i];  
                                                 // from KK model in flux form;
      S2[i] = (rho[i]*vw - Q[i]) / kkl.tau; // all interaction in relaxation
    }

  } // end KK


  /***************************************************************/
  /** Macroscopic version of the IDM                             */
  /***************************************************************/

  else if (choice_model==3)  
  {

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<TINY_VALUE) rho[i]=TINY_VALUE;
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]=rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }

    // calculate advanced rho and v fields (antic_factor=1 in Mac. IDM!!)
    // with  shift_in_x_mca ( -> xa=x-1/rho) instead of shift_in_x!

    shift_in_x_mca (rho, rho,  rhodelta);
    shift_in_x_mca (rho, v,    vdelta);

    for (i=0; i<=nx; i++)
    {
      double vdiff     = v[i] - vdelta[i];
      double xdiff     = 0.5*(1/rho[i]+1/rhodelta[i]) - 1/rhomax;
      double dxstarmin = idm.s0;
      double dxstar    = idm.s0 + Tr_loc[i] * v[i] 
                          + 0.5 * v[i]*vdiff * sqrt(tau0/(idm.b*v0_loc[i]));
      if(dxstar<dxstarmin) dxstar=dxstarmin;
      double a_IDM     = v0/tau0 
            * (1 - pow(v[i]/v0_loc[i], idm.delta) - SQR(dxstar/xdiff) );

      // output of calc_rhs

      F1[i] = Q[i];
      F2[i] = rho[i] * SQR(v[i]);               // equiv. to mik
      S2[i] = rho[i]* a_IDM;
    }   

  } // end Mac. IDM model

  // *********************************************************************
  //                  Macroscopic version of the OVM-Model      
  // ********************************************************************


  else if (choice_model==4)
  {
    double deltave_mca;         // anticipated equilibrium velocity for 
    double rho_shift[NXMAX+1];   // average of rho, rhodelta

    for (i=0; i<=nx; i++)
    {

    // operations of "test_and_calc"

      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]=rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];

    }

    // calculate advanced rho and v fields
 
    shift_in_x_mca (rho, rho,  rhodelta);

    for (i=0; i<=nx; i++) rho_shift[i]=.5*(rho[i]+rhodelta[i]); 
    shift_in_x_mca (rho_shift ,rho,  rhodelta);
    for (i=0; i<=nx; i++)
    {
      // ve(rho_a) with ve(rho)=ve_Bando(rho) or max(v0,s/T) => equil.cc
      deltave_mca = intp(veqtab, NRHO+1, rhodelta[i],0,rhomax); 
      F1[i] = Q[i];
      F2[i] = rho[i] *SQR(v[i]);
  
      S2[i] = rho[i]*(1/tau0)*(deltave_mca-v[i]) ;
      if (int(S2[i]) < 0) S2[i]=0.0;   
    }


 

  }
//End OVM





//####################### Black scholes ###################

/*
    Linear Black-Scholes PDE: linear PDE with space-dependent coefficients;
    
     ----------------------------------------------------
     | d_t C = -r0 C + r0 x d_x C + 1/2 sig^2 x^2 d_x^2 C |
     ----------------------------------------------------

    rho here: Price C of option in 1000 DM (since rho mult. by 1000 on outp)
    x   here: Stock price at time t0 in 1/1000 DM since divided by 1000 on outp
    t   here: time t0 (< strike time t_s)  (PDE backwards in time)
        time in years/60, so that output is in years (usually s vs. min)

   Notice: scaling of eqs only necessary for t: 
   r with 1/t, sig with 1/sqrt(t)
   No scaling of rho(since linear), or x (x and d_x always in same power)

   In conservative form we rewrite the BS eq.  as

     ----------------------------------------------------------------
     | d_t C = (sig^2-2 r0) C + d_x[ (r0-sig^2)xC + 1/2 sig^2 x^2d_x C] |
     ----------------------------------------------------------------
->   F1 = (sig^2-r0) x C - 1/2 sig^2 x^2d_x C
     S1 = (sig^2-2 r0) C
*/

  else if (choice_model==5)
  {
    if(choice_method==1){
      error(" Black scholes is  parabolic; use McCormack instead of upwind");
    }

    double sig = v0;          // "volatility" = sqrt(variance) in 1 year
    double r0   = Tr;                //  risk-free interest rate (per year)

    double xmin = tau0;               // stock price at left boundary xmin
                                 // max stock price = xmin+xmax (of .inp file)

    int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                        // ishift=1 for downwind diffrences



    for (i=0; i<=nx; i++)
      if (rho[i]<TINY_VALUE) rho[i] = TINY_VALUE;
    //    if (rho[i]<TINY_VALUE) error("calc_rhs for Black-Scholes: Prize rho<0");


    for (i=0; i<=nx; i++)
    {
      double x = xmin + xmax*i/nx; // induces nonlinearity of BSE !!!

      F1[i] = (SQR(sig)-r0) * x * rho[i] 
	  - 0.5 * SQR(sig*x) * ( rho[i+ishift]-rho[i+ishift-1])/dx;
      S1[i] = (SQR(sig)-2.*r0) * rho[i];  //!!

      F2[i] = 0;
      S2[i] = 0.;  
    }



  }//End Black-Scholes

  //################# 3 Phase macro model! (mar08)
  else if (choice_model==8){
    mac3phases.calc_rhs(rho,Q,F1,F2,S1,S2,  choice_BC, show_calc_rhs);
  }

//####################### Fokker Planck equation #############

  // FPE: dP(x,t)/dt = - d/dx (A P) + 1/2 d^2/dx^2 (B P)
  // P(x,t)===rho(x,t)

  else if (choice_model==50){

    fpe.calc_rhs(rho,F1,S1,downwind_diff);

    /*
    if(choice_method==1){
      error(" FPE is parabolic; use McCormack instead of upwind");
    }

    double fpeA=Tr;
    double fpeB=v0;   // B = 2D (B=Kramers-M coeff, D=diffusion coefficient)

    int    ishift = downwind_diff ? 1 : 0; // Only for flow derivatives:
                                        // ishift=1 for downwind diffrences


    for (i=0; i<=nx; i++){
      F1[i]= fpeA* rho[i]
              - 0.5 * fpeB * ( rho[i+ishift]-rho[i+ishift-1])/dx;
      S1[i]=0;

      F2[i] = 0;
      S2[i] = 0.;  
    }
    
    */


  }
  //End FPE

//################ InnovationFokker Planck equation ###################

  // Inno-FPE: dP(x,t)/dt = lambda*P(x,t) * (x-erw(x)) + d^2/dx^2 (D P)
  // P(x,t)===rho(x,t)

  else if (choice_model==51){
    fpe_innov.calc_rhs(rho,F1,S1,downwind_diff);
  }
  //End FPE_innov

  else if (choice_model==21){

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<TINY_VALUE)  rho[i] = TINY_VALUE;
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]=rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }

    vmm.calc_rhs(downwind_diff,rho,Q,F1,F2,S1,S2);
  }
  //End VMM

  else if (choice_model==10){

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<TINY_VALUE)  rho[i] = TINY_VALUE;
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]=rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }
    sgm.calc_rhs(downwind_diff,rho,Q,F1,F2,S1,S2);
  }
  //End SGM


  else error(" Sorry, model for given choice-model not implemented");

  //################################################################
  // test output
  //################################################################


  if (show_calc_rhs){
    int ix_show_min=0;
    int ix_show_max=5;
    //int ix_show_max=nx;
    cout <<"\ncalc_rhs: downwind_diff="<<downwind_diff
	 <<" debug output:"<<endl;
    printf("ix\t  rho\t Q\t V\t F2    \t S2  \n");
    for (i= ix_show_min; i<= ix_show_max; i++){
        printf("%i\t %.4f\t %.4f\t %.2f\t %.2f\t %.4f\n",
	       i, rho[i], Q[i], Q[i]/rho[i], 
                 F2[i],     
                 S2[i]);
    }
  }



  //################################################################
  // calculate sources from on- and off-ramps for all models above 
  //################################################################

    // only Lee-KK with Gaussian ramp flow density!!
  int gauss = (choice_model==2) ? true : false; 


  if((choice_model!=5)&&(choice_model<50)){
    for(i=0;i<=nx; i++) S1[i]=0;
  }

  if(choice_rmp==1)
  {
    add_rmp (S1, S2, rho,Q, x_rmp,dx_rmp,flow_rmp, gauss );
         // flow_rmp interpolated from the array Q_rmp[] in timestep 
  }

  if((choice_rmp==2)||(choice_rmp==3))
  {
    int i_rmp;
    for (i_rmp=0; i_rmp<=n_rmps; i_rmp++)
      add_rmp (S1, S2, rho, Q, x_rmps[i_rmp],dx_rmps[i_rmp],
               flow_rmps[i_rmp], gauss );
         // flow_rmps[] interpolated from the array Q_rmps[][] in timestep 
    //cout <<"add_rmp: i_rmp="<<i_rmp<<" flow_rmps[i_rmp]="<<flow_rmps[i_rmp]<<endl;
  }
  if(choice_rmp==3)
    {
    int i_rmp;
    for (i_rmp=0; i_rmp<=n_rmps_c; i_rmp++)
      add_rmp (S1, S2, rho, Q, x_rmps_c[i_rmp],dx_rmps_c[i_rmp],
               flow_rmps_c[i_rmp], gauss );  
    }
  
 
}

/*********************************************************************/
/****** end calc_rhs                                              ****/
/*********************************************************************/




/*********************************************************************/
/*  function add_rmp                                                */
/********************************************************************/

void add_rmp (double S1[], double S2[], double rho[], double Q[],
              double x_rmp, double dx_rmp, double flow_rmp, int gauss)

  /* adds the effects of one ramp to the source terms S1 and S2
     The ramp is at position x_rmp with merging length dx_rmp; its flow
     flow_rmp is >0 for an on-ramp, <0 for an off-ramp.
     The velocity of the merging or leaving cars is the local velocity
     of the main road.    
     If gauss is on, the ramp flow distribution is Gaussian centered at x_rmp
     with sqrt(variance)=L; otherwise it is constant in [xrmp-L/2,xrmp+L/2]
     */

    { 
      if (!gauss)    // constant ramp flow density
      {
        int    nx_rmp  = (int)(dx_rmp/dx);  // (actual length = nx_rmp*dx)
        int    ix_min  = (int)((x_rmp-0.5*dx_rmp)/dx);
        int    ix_max  = ix_min + nx_rmp; 
        int    ix,nx   =(int)(xmax/dx);
   
        if( ix_min<1)  ix_min=1;
        if( ix_max>nx-1) ix_max=nx-1;

        for(ix=ix_min+1;ix<=ix_max; ix++) 
        {
          double source = flow_rmp/(nx_rmp*dx);
          S1[ix] += source;
          S2[ix] += source*Q[ix]/rho[ix];    // merging with velocity of main road 
        }
      }
      else        // Gaussian ramp density
      {
        int    nx_rmp  = (int)(6.*dx_rmp/dx);          //  +- 3 sigma
        int    ix_min  = (int)((x_rmp-3*dx_rmp)/dx);
        int    ix_max  = ix_min + nx_rmp; 
        int    ix,nx   =(int)(xmax/dx);

        if( ix_min<1)  ix_min=1;
        if( ix_max>nx-1) ix_max=nx-1;

        double prefac = 1./( dx_rmp*sqrt(2*PI));
        for(ix=ix_min+1;ix<=ix_max; ix++) 
        {
          double xdev = ix*dx - x_rmp;
          double gauss_dist = prefac*exp(-0.5*SQR(xdev/dx_rmp));
          double source = flow_rmp * gauss_dist;
          S1[ix] += source;
          S2[ix] += source*Q[ix]/rho[ix];    // merging with velocity of main road 
	  //cout <<"prefac="<<prefac<<" exp="<<exp(-0.5*SQR(xdev/dx_rmp))
          //     <<" x="<<ix*dx<<" source="<<source<<endl;
        }

      }
    }



double  cont_rmpflow(double rmpflow_in,double flow_main,double flow_crit, 
		     double dens_main,double *stor_ptr)

// Implements a simple on-ramp-flow control by cutting off the total flow
// = flow_main + rmpflow_in at some critical value flow_crit
// The additional parameter dens_main (density) can be used to devise
// more elaborated controls. flow_main and dens_main should be measured
// 100m..500m upstream of the controlled on-ramp.
// This implementation does not prohibit a negative effective ramp flow
// corresponding to a controlled deviation at a neighbouring off-ramp.
// Macroscopically, the whole thing can be considered 
// as simulation of a controlled intersection.

{

  double calc_flow=0;       // allowed ramp outflow given by on-ramp control
  double Qout=0.15;         // outflow (veh/s/lane) of queuing vehicles
  double queued_veh = *stor_ptr; // how many vehicles  are waiting?  
  double flow_wanted = (queued_veh<=0) ? rmpflow_in : Qout;

  calc_flow = (flow_crit-flow_main)<flow_wanted
            ? flow_crit - flow_main
            : flow_wanted;

  *stor_ptr = *stor_ptr + (rmpflow_in - calc_flow)*dt;

  //printf("%lf\t%lf\n",*stor_ptr,flow_main);
  
  return(calc_flow);
}



void dyn_v0(double rho[])

{
double v_tab[NBC];
int i;

for(i=0; i<=n_jmps_v0; i++)
	{
	  double rho_loc = intp(rho, int(xmax/dx)+1,x_v0[i] ,0,xmax);
	  
	  if (rho_loc > thrsh_x[i]) v_tab[i] = (v0_x[i] - delta_v0_x[i]);
	  else v_tab[i] = v0_x[i];
	}
for (i=0; i<=(int)(xmax/dx); i++)
  {
    //!! Factor trucks.v0TruckCorr here relevant???
       v0_loc[i] = intpextp(x_v0,v_tab,n_jmps_v0+1,i*dx);
       //printf("%d\t%lf\n",i,v0_loc[i]);
       
  }
}



