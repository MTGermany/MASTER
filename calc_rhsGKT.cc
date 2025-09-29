

// WATCH OUT: calc_rhsGKT taken only
// if testIsolatedGKT==true in calc_rhs
// but also calc_rhsGKT uses normal helpers add_rmp, shift_in_x etc
// furthermore, add_rmp is
// ALWAYS added at the end of calc_rhs in calc_rhs.cc


//###########################################
//  function calc_rhsGKT                  
//###########################################


void calc_rhsGKT (int choice_model, bool downwind_diff,
              double rho[],double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[], 
              bool debug)

  /* Calculates the right-hand sides (fluxes Fi and sources Si)
     of macroscopic traffic models of the form
              ________________________________________________
              |  d_t rho(x,t)  = -d_x F1 + S1 (ramp flows)    |
              |  d_t Q(x,t)    = -d_x F2 + S2                 |
              ------------------------------------------------

     Input: 
       - choice_model: defines the traffic model; choice_model==0: GKT model
       - downwind_diff: whether downwind differences should be used
         for calculating derivatives (only relevant if fluxes with gradients)
       - The fields rho,Q. They will be modified only, if 
         something unphysical happens
       - switch "debug" to show some debugging information

     Output: 
       - The right-hand sides of the PDE:
         fluxes F1,F2, and sources S1,S2 of the model.
       - The sources include ramps: S1 = ramp flow/merging length 
                                    S2+= S1*(velocity when merging) 
       - Flow-conserving bottlenecks realized by gradients of the
         parameters T (=Tr) or v0, given by the (global) arrays
         Tr_loc and v0_loc

     Differences from the inline function calc_rhs in calc_rhs.cc 
     (which is used if calc_rhs.testIsolatedGKT=false): (MT 2025)
     All changes marked by MT 20xx

     Tested: Below rhomax/2 no visible differences, 
     particularly not in the istability region (MT 2025-09)

       - [the density is not only bounded from below (rho>0) 
         but also from above (rho<rhomax)]

       - [shift_in_x(.) is the same as shift_in_xGKT]

       - (1) change of the relaxation time tau0->tau1 [book ed. 2, (10.30)]
         to avoid high-dens local instabilities

       - (2) Some small high-density diffusion directly in the conservative
         terms F1 and F2, D=Dmax*(1-ve/v0)^2, Dmax=10 m^2/s (very small)

 
   */


{

  //cout <<" in calc_rhsGKT, choice_model="<<choice_model<<endl;
  
  if(false){
    cout<<" in calc_rhsGKT: downwind_diff="<<((downwind_diff) ? "true" : "false")<<endl;
  }

  int    i;
  int    nx   =(int)(xmax/dx);
  double v[NXMAX+1];                         // velocity
  double rhodelta[NXMAX+1];                  // nonlocal density
  double vdelta[NXMAX+1];                    // nonlocal velocity

  
  // #########################################################
  // Standard GKT model with new reliable highspeed corr tau0->tau>tau0
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
  // ACHTUNG only if testIsolatedGKT==true in calc_rhs; otherwise
  // code in calc_rhs taken !!!
  // #############################################################

  if  ( choice_model!=0){
    cerr<<" Error: Only the GKT (choice_model=0) contained in this version"
	<<endl;
    exit(-1);
  }
  if(choice_model==0){


    double Dplusmax=10;  // for later high-dens diffusion MT 2025
    const double durch_dx = 1./dx; // for later high-dens diffusion MT 2025
    
    // Check input (rho>0, Q>=0) and calculate V

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<TINY_VALUE)        rho[i] = TINY_VALUE;
      if (rho[i]>rhomax-TINY_VALUE) rho[i] = rhomax-TINY_VALUE; //MT 2016
      if (Q[i]<rho[i]*TINY_VALUE)   Q[i]   =rho[i]*TINY_VALUE;
      v[i] = Q[i]/rho[i];
    }

    // calculate advanced rho and v fields rhodelta[], vdelta[] 
    
    shift_in_x (antic_factor, v, rho, rhodelta);
    shift_in_x (antic_factor, v, v,   vdelta);

    // calculate rhs of flow equation F1, F2 and S2
    // S1, S2inh calculated in ramp-source section below

    for (i=0; i<=nx; i++)
    {
      double sqrtA     = intp(sqrtAtab, NRHO+1, rho[i], 0, rhomax);
      double denom     = 1-rhodelta[i]/rhomax; // rhodelta<rhomax-TINY_VALUE
      double dvrel     = (v[i]-vdelta[i]) / (sqrtA * sqrt2 * v[i]  );   
      double bolzm_fact= ( (dvrel>DVMIN) && (dvrel<DVMAX) )
                       ? intp( Btab, NDV+1, dvrel, DVMIN, DVMAX)
                       : (dvrel>0) ? 2.*SQR(dvrel) : 0.; 



      // calculate F1, F2 and S2; S1 calculated in ramp-source section
      // below

      F1[i] = Q[i];
      F2[i] = rho[i] * SQR(v[i]) * (1.+ SQR(sqrtA)); // hat wenig Einfluss

      // change of the relaxation time to avoid high-dens local instabilities 
      // MT 2025
      
      double ve=max(0.0001*v0_loc[i], intp(veqtab, NRHO+1, rho[i], 0, rhomax));
      //double tau1=tau0;
      double tau1=max(tau0,2*dt*(1+2*v0_loc[i]/ve)); //MT 2025



     

      double S2free = (rho[i]*v0_loc[i] - Q[i]) / tau1; // MT 2025 tau0->tau1
      double S2brake=  rho[i]*v0_loc[i] // rhodelta better but not consistent
	  * SQR( Tr_loc[i]*sqrtA*rhodelta[i]*v[i]) * bolzm_fact
	  / (SQR(denom)*tau1*Arhomax);  // MT 2025 tau0->tau1
      S2[i] = S2free  - S2brake;


    // include high-density diffusion directly onto F1 and F2 (MT 2025)

      double D_corr      = Dplusmax * SQR(1-ve/v0);      // additl diffusion

      if (downwind_diff){
	 F1[i] += - D_corr * ( rho[i+1] - rho[i]) * durch_dx;
         F2[i] += - D_corr * ( Q[i+1]   - Q[i]  ) * durch_dx;
      }

      else{
	 F1[i] += - D_corr * ( rho[i] - rho[i-1]) * durch_dx;
         F2[i] += - D_corr * ( Q[i]   - Q[i-1]  ) * durch_dx;
      } 
    


      if(false){
	//if(i==nx/2){
	cout<<" rho[i]="<<rho[i]<<" v[i]="<<v[i]<<" ve="<<ve<<" denom="<<denom<<" tau1="<<tau1<<" tau1/tau0="<<tau1/tau0
		      <<" downwind_diff="<<downwind_diff
		      <<" D_corr="<<D_corr<<endl;
      }
      
    }

 

    // test output
    if (debug)
    {
      printf("\ncalc_rhs: Output \n\n");
      printf("ix\t  rho\t F2=E     \t v0_loc   \t S2  \n");
      for (i= ix_show_min; i<= ix_show_max; i++)
      {
        printf("%i\t %.6f\t %.6f\t %.6f\t %.6f\n",
                 i,  
                 rho[i],     
                 F2[i],     
                 v0_loc[i], 
                 S2[i]);
      }
    }

  }// end GKT calculation with new high-dens corr



}

//########################################
// end calc_rhsGKT
//########################################




