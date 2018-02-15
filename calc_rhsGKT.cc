


//########################################
//  function add_rmpGKT                                                */
//########################################


void add_rmpGKT (double S1[], double S2[], double rho[], double Q[],
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
          S2[ix] += source*Q[ix]/rho[ix];  
        }

      }
}


//#################################################################
//  function shift_in_xGKT  
//#################################################################


void shift_in_xGKT (double antic_factor, 
     const double v[], const double f[], double fshifted[])

  //   Calculates, from the input array "f" representing f(x),
  //   an array "fshifted" representing f(x+dx_shift) 
  //   The shifted distance is s=antic_factor * (1/rhomax + Tr*v). 
  //   The extrapolation to the right depends on the BC:
  // !! changes here imply changes in function "boundary_vals" and vice versa

{
  double dx_shift,rest;
  int    i, idx_shift, nx=(int)(xmax/dx);

  if(antic_factor > 0)
  {
    for(i=0; i<=nx;i++)
    {
      double T_xt = Tr_loc[i]*trucks.TrTruckCorr;
      dx_shift  = antic_factor * ( 1/rhomax + T_xt *v[i] );
      idx_shift = (int)(dx_shift/dx);
      rest   = dx_shift/dx - idx_shift;

      if((i+idx_shift+1)<=nx)
        fshifted[i]    = (1-rest)*f[i+idx_shift] + rest*f[i+idx_shift+1];

      // extrapolation depending on the BC: 0=periodic,1=free,2=Neumann, 3=Dirichlet

      else 
      { 
        if      (choice_BC==0) fshifted[i] 
            = (1-rest)*f[i+idx_shift-nx+1] + rest*f[i+idx_shift-nx+2];
        else if ((choice_BC==1) || (choice_BC==4)) fshifted[i] 
            = f[nx] + (dx_shift/dx - nx+i) * (f[nx]-f[nx-1]);
        else fshifted[i] 
	    = f[nx];  // const extrapolation in case of Neumann  or fixed BC

      }
    }
  }
  else for(i=0; i<=nx;i++) fshifted[i]=f[i];
}

//#################################################################
//  function shift_in_xLWR  
//#################################################################


void shift_in_xLWR (double antic_factor, 
     const double f[], double fshifted[])

  //   Calculates, from the input array "f" representing f(x),
  //   an array "fshifted" representing f(x+dx_shift) 
  //   The shifted distance is ds=antic_factor * grid size dx 
  //   The extrapolation to the right depends on the BC:
  // !! new BC in function "boundary_vals" may imply changes here!

{
  int nx=(int)(xmax/dx);

  if(antic_factor > 0){
    for(int i=0; i<=nx;i++){
      double dx_shift  = antic_factor * dx;
      int idx_shift = (int)(dx_shift/dx);
      double rest   = dx_shift/dx - idx_shift;

      if((i+idx_shift+1)<=nx)
        fshifted[i]    = (1-rest)*f[i+idx_shift] + rest*f[i+idx_shift+1];

      // extrapolation depending on the BC: 0=periodic,1=free,2=Neumann, 3=Dirichlet

      else{ 
        if      (choice_BC==0) fshifted[i] 
            = (1-rest)*f[i+idx_shift-nx+1] + rest*f[i+idx_shift-nx+2];
        else if ((choice_BC==1) || (choice_BC==4)) fshifted[i] 
            = f[nx] + (dx_shift/dx - nx+i) * (f[nx]-f[nx-1]);
        else fshifted[i] 
	    = f[nx];  // const extrapolation in case of Neumann  or fixed BC

      }
    }
  }
  else for(int i=0; i<=nx;i++) fshifted[i]=f[i];
}



//#######################################################
//  function  intpGKT  
//#######################################################

double intpGKT (const double tab[], int n,double x, 
	     double xmin, double xmax)
   /* intp interpolates the array tab with n 
      equidistant points
      in [xmin,xmax] at the location x; an error message is produced on
      attempting extrapolation */
{
  int nloc=n-1; 
  double intp_value;
  double ir   = nloc*(x-xmin)/(xmax-xmin);
  int    i    = (int) ir;
  double rest = ir-i;
  if ((i>=0) && (i<=nloc-1))  intp_value =  (1-rest) * tab[i] + rest*tab[i+1];
  else if (i==nloc) intp_value = tab[nloc];
  else {
     cout << "intp: index i = "<<i<<" (ir="<<ir<<") out of range\n";
     intp_value = tab[nloc];
  }

  return(intp_value);
}


//###########################################
//  function calc_rhsGKT                  
//###########################################


void calc_rhsGKT (int choice_model, bool downwind_diff,
              double rho[],double Q[],  
              double F1[], double F2[],  
              double S1[], double S2[], 
              int debug)

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

     Differences from the inline function in calc_rhs.cc 
     (which is used if calc_rhs.testIsolatedGKT=false): 

       - the density is not only bounded from below (rho>0) 
         but also from above (rho<rhomax) (marked as all the rest by //MT 2016)
       - counter vanishing relax. times [book, (9.45)] by multiplying the 
         flow relaxation term S2 by (1-rho/rhomax)
 
   */


{

  if(false){
    cout<<" in calc_rhsGKT: downwind_diff="<<((downwind_diff) ? "true" : "false")<<endl;
  }

  int    i;
  int    nx   =(int)(xmax/dx);
  double v[NXMAX+1];                         // velocity
  double rhodelta[NXMAX+1];                  // nonlocal density
  double vdelta[NXMAX+1];                    // nonlocal velocity

  double tiny_value=1.e-6;
  double rangeBoltzmann=3;
  
  // #########################################################
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
  // ACHTUNG Berechnung bei calc_rhs, nicht diese,  wird genommen !!!!
  // #############################################################

  if  ( (choice_model!=0) && (choice_model!=9)){
    cerr<<" Error: Only the GKT (choice_model=0) and a LWR mockup (choice_model=9)"
	<<endl<<" contained in this version"<<endl;
    exit(-1);
  }
  if(choice_model==0){

    // Check input (rho>0, Q>=0) and calculate V

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<tiny_value)        rho[i] = tiny_value;
      if (rho[i]>rhomax-tiny_value) rho[i] = rhomax-tiny_value; //MT 2016
      if (Q[i]<rho[i]*tiny_value)   Q[i]   =rho[i]*tiny_value;
      v[i] = Q[i]/rho[i];
    }

    // calculate advanced rho and v fields rhodelta[], vdelta[] 
    
    shift_in_xGKT (antic_factor, v, rho, rhodelta);
    shift_in_xGKT (antic_factor, v, v,   vdelta);

    // calculate rhs of flow equation F1, F2 and S2
    // S1, S2inh calculated in ramp-source section below

    for (i=0; i<=nx; i++)
    {
      double sqrtA     = intp(sqrtAtab, NRHO+1, rho[i], 0, rhomax);
      double denom     = 1-rhodelta[i]/rhomax;
      double dvrel     = (v[i]-vdelta[i]) / (sqrtA * sqrt2 * v[i]  );   
      //double dvrel     = (v[i]-vdelta[i]) / (sqrtA * v[i]  );   
      //double dvrel     = 0;
      double bolzm_fact= ( (dvrel>-rangeBoltzmann) && (dvrel<rangeBoltzmann) )
                       ? intp( Btab, NDV+1, dvrel, -rangeBoltzmann, rangeBoltzmann)
                       : (dvrel>0) ? 2.*SQR(dvrel) : 0.; 
      // calculate F1, F2 and S2; S1 calculated in ramp-source section
      // below

      F1[i] = Q[i];
      F2[i] = rho[i] * v[i]*v[i] * (1.+ sqrtA*sqrtA);
      //F2[i] = rho[i] * SQR(v[i]) * (1.);   // !!!
      double S2free = (rho[i]*v0_loc[i] - Q[i]) / tau0;
      double S2brake= v0_loc[i] * rhodelta[i]
                      * SQR( Tr_loc[i]*sqrtA*Q[i]) * bolzm_fact
	              / (SQR(denom)*tau0*Arhomax);
      S2[i] = S2free  - S2brake;
      S2[i] *= denom;  //MT 2016 counter vanishing relax. times (9.45)
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

  }

  else if(choice_model==9){ // LWR mockup with parabolic FD (MT 2016)
    // flux term F2 of flow-conservative formulation = Q^2/rho (book (9.30))
    // source term S2 of flow-conservative formulation = (rho*Ve-Q)/tau0 + S2inh


    // Check input (rho>0, Q>=0) and calculate V

    for (i=0; i<=nx; i++)
    {
      if (rho[i]<tiny_value)        rho[i] = tiny_value;
      if (rho[i]>rhomax-tiny_value) rho[i] = rhomax-tiny_value;
      if (Q[i]<rho[i]*tiny_value)   Q[i]   =rho[i]*tiny_value;
      v[i] = Q[i]/rho[i];
    }

    // calculate advanced v field (just to realize downwards dependence for num stability)
    
     shift_in_xLWR (antic_factor, rho,   rhodelta);

    // calculate rhs of flow equation F1, F2 and S2
    // S1, S2inh calculated in ramp-source section below

    for (i=0; i<=nx; i++)
    {
      double ve=intp(veqtab, NRHO+1, rhodelta[i], 0, rhomax);
      //double ve=v0_loc[i]*(1-rhodelta[i]/rhomax);
      ve=min(v0_loc[i], max(0., ve));

      F1[i] = Q[i];
      F2[i] = rho[i] * SQR(v[i]);
      S2[i] = (rho[i]*ve-Q[i])/tau0;
    }

  }


  // end actual GKT or LWR calculation (choice_model=0 or =9)



  //###########################################
  // calculate sources from on- and off-ramps for all models 
  //###########################################

 // initialize
  for(i=0;i<=nx; i++){
    S1[i]=0;
  }

      // only Lee-KK with Gaussian ramp flow density!!
  int gauss = (choice_model==2);

 
  if(choice_rmp==1)
  {
    add_rmpGKT (S1, S2, rho,Q, x_rmp,dx_rmp,flow_rmp, gauss );
         // flow_rmp interpolated from the array Q_rmp[] in timestep 
  }

  if((choice_rmp==2)||(choice_rmp==3))
  {
    int i_rmp;
    for (i_rmp=0; i_rmp<=n_rmps; i_rmp++)
      add_rmpGKT (S1, S2, rho, Q, x_rmps[i_rmp],dx_rmps[i_rmp],
               flow_rmps[i_rmp], gauss );
         // flow_rmps[] interpolated from the array Q_rmps[][] in timestep 
  }
  if(choice_rmp==3)
    {
    int i_rmp;
    for (i_rmp=0; i_rmp<=n_rmps_c; i_rmp++)
      add_rmpGKT (S1, S2, rho, Q, x_rmps_c[i_rmp],dx_rmps_c[i_rmp],
               flow_rmps_c[i_rmp], gauss );  
    }
  
 
}

//########################################
// end calc_rhsGKT
//########################################




