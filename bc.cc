
/********************************************************************/
/*  function  boundary_vals                                        */
/********************************************************************/
 
void   boundary_vals(double field[], int nx, 
                     int choice_BC, double left, double right)

// Defines the two boundary points field[0], field[nx] of "field" 
// !! if new BC are defined, the function "shift_in_x","shift_in_x_mca"
//  must be modified accordingly and so must "timestep", "input".

// flow notice: whether steady-state flow or flow from files is assumed
// is determined in timestep.cc possibly setting Qright, Qleft to the
// steady-state values
 
{
  if (choice_BC==0)      /* periodic BC */
  {
    field[0]  = field[nx-1];
    field[nx] = field[1];
  }

  else if(choice_BC==1)  // Dirichlet (fixed) BC at inflow, 
                         // lin. extrapolation ("free BC") at outflow
  { 
    field[0]  = left; 
    field[nx] = 2.*field[nx-1] - field[nx-2];
  }

  else if(choice_BC==2)   // Dirichlet (fixed) BC for rho at inflow, Q equil.
                          // homogeeous Neumann BC to the  right
  {
    field[0]  = left;
    field[nx] = field[nx-1];
  }

  else if(choice_BC==3)   // Dirichlet BC both at inflow and outflow,
                          // Q in equil. (this is selected in timestep.cc)
  {
    field[0]  = left;
    field[nx] = right;
  }

  else if(choice_BC==4)   // Free BC both at inflow and outflow
  {
    field[0]  = 2.*field[1]    - field[2];
    field[nx] = 2.*field[nx-1] - field[nx-2];
  }
  else if(choice_BC==5)  // Dirichlet BC for rho,Q at inflow and outflow
  {
    field[0]  = left;
    field[nx] = right;
  }
  else if(choice_BC==6)  // homogeneous Neumann BC to the left and right
  {
    field[0] = field[1];
    field[nx] = field[nx-1];
  }

  else if(choice_BC==7)   // Dirichlet (fixed) BC for rho,Q at inflow, 
                          // homogeneous Neumann BC to the  right
  {
    field[0]  = left;
    field[nx] = field[nx-1];
  }

  else{;}   // choice_BC=8,9,... handled by following boundary_vals() function

    //error("boundary_vals: choice_BC>7 or <0 not implemented here");
}

/********************************************************************/
/*  traffic-dynamic specific function  boundary_vals                */
/********************************************************************/

void   boundary_vals(double rho[], double Q[], int nx, 
                     int choice_BC, double rholeft, double Qleft,
                     double rhoright, double Qright)

  // Defines the four boundary points rho[0], rho[nx], Q[0], Q[nx]
  // implements specific dynamic BC of traffic-dynamic equations
  // (choice_BC=8,9,10) 
  // !! if new BC are defined, the function "shift_in_x","shift_in_x_mca"
  //  must be modified accordingly and so must "timestep", "input"

{

  // check input
  const double epsilon=0.0001;

  if(rholeft<epsilon)          rholeft=epsilon;
  if(rholeft> rhomax-epsilon)  rholeft=rhomax-epsilon;
  if(rhoright<epsilon)         rhoright=epsilon;
  if(rhoright> rhomax-epsilon) rhoright=rhomax-epsilon;

  if(Qleft <v0*epsilon)        Qleft=v0*epsilon;
  if(Qright<v0*epsilon)        Qright=v0*epsilon;

  // simple boundary conditions

  if((choice_BC>=0) && (choice_BC<=7) )
  {
    boundary_vals(rho,nx,choice_BC,rholeft,rhoright);
    boundary_vals(Q,  nx,choice_BC,Qleft,  Qright);
  }

  // dynamic BC treated here; actual distinction: timestep.cc!!
  // choice_BC  upstream BC                  downstream BC
  // -------------------------------------------------------------
  // 8          Dir_rho,  Dir_Q              Dir_rho, equil_Q
  // 9          equilFree_rho(Q), Dir_Q      equilCong_rho(Q), Dir_Q
  // 10         Dir_rho,  Dir_Q              equil_rho(v), equil_Q(v)

  // Distinctions between 8-10 in timestep.cc

  else if( (choice_BC==8)  // dynamic BC using rho and Q from .BC file
	|| (choice_BC==9) // dynamic BC using Q from.BC file,
                           // and rho=inverse of Qe(rho) -> timestep
        || (choice_BC==10)) // dynamic BC as 8, but using downstream velocity
                            // -> rho_down=inverse of Ve(rho), Q=Qe(rho)
  {

    // Determine whether Dirichlet or hom. Von-Neumann BC are used

    const double alpha   = 0.95;  // fit factor (e.g., 0.95) 
    const double Qreduce = 0.98;  // fit factor (e.g., 0.98)
    bool free_traffic    = ((rho[0]<alpha*rhoQmax) && (rho[2] <alpha*rhoQmax));
    bool jam_dissolve    = (Qreduce*Q[2]-Qleft  > 0);
    bool dirichlet_up    = free_traffic  || jam_dissolve;
    bool dirichlet_down  = (choice_method==1) 
        ? true                                           // upwind scheme
        :( (Q[nx-1]-Qright) * (rho[nx-1]-rhoright) < 0); // vg<0; for McCormack

    // ######### Apply corresponding BC

    if (dirichlet_up){
      rho[0]  = rholeft;
      Q[0]    = Qleft;
    }

    else{
      //if (it%50 == 0) cout<<"N!! ";  // Quick-Docu hack
       rho[0]  = rho[1];
       Q[0]    = Q[1];
    };

    if (dirichlet_down){
      rho[nx]  = rhoright;
      Q[nx]    = Qright;
    }
    
    else{
      rho[nx]  = rho[nx-1];
      Q[nx]    = Q[nx-1];
    };

  }  // end choice_BC==8-10

  else  error("boundary_vals: choice_BC>10 or <0 not implemented");

  if( show )
    cout << "boundary_vals: " << endl
           <<" rholeft="   <<rholeft  
           <<" \t\t Qleft="  <<Qleft    
           <<" \t rhoQmax="<<rhoQmax << endl
	   <<" rho[0]="    <<rho[0]   
	   <<" \t rho[1]=" <<rho[1]   
	   <<" \t Q[0]="   <<Q[0]     
	   <<" \t Q[1]="   <<Q[1]   <<endl   
           <<" rhoright="  <<rhoright
           <<" \t\t Qright=" <<Qright  << endl
           <<" rho[nx]="   <<rho[nx] 
           <<" \t Q[nx]="  <<Q[nx]   
           <<" \t\t rho[nx-1]="<<rho[nx-1] 
           <<" Q[nx-1]="<<Q[nx-1] 
           <<endl;

}





//####################################################################
//  function shift_in_x  
//####################################################################


void shift_in_x (double antic_factor, 
     const double v[], const double f[], double fshifted[])

  //   Calculates, from the input array "f" representing f(x),
  //   an array "fshifted" representing f(x+dx_shift) 
  //   The shifted distance is s=antic_factor * (1/rhomax + Tr*v). 
  //   The extrapolation to the right depends on the BC:
  // !! changes here imply changes in functions
  // !! "shift_in_x_mca" and "boundary_vals" and vice versa

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

      // extrapolation depending on the BC: periodic,free,Neumann=Dirichlet

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



//####################################################################
//  function shift_in_x_mca  
//####################################################################



void shift_in_x_mca (const double rho_t[], 
                     const double f[], double fshifted[])

  //   Calculates, from the input array "f" representing f(x),
  //   an array "fshifted" representing f(x+dx_shift) 
  //   The shifted distance is s=(1/rho) 
  //   The extrapolation to the right depends on the BC:
  // !! changes here imply changes in function "boundary_vals" 
  //  and vice versa

{
  int nx=(int)(xmax/dx);

  if(antic_factor > 0){
    for(int i=0; i<=nx;i++){
      double dx_shift= 1./rho_t[i];      // only this differs from shift_in_x
      //double dx_shift  = antic_factor * dx;
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

//#################################################################
//  function shift_in_xLWR  (only choice_model==9)
//#################################################################


void shift_in_xLWR (double antic_factor, 
     const double f[], double fshifted[])

  //   Calculates, from the input array "f" representing f(x),
  //   an array "fshifted" representing f(x+dx_shift) 
  //   The shifted distance is ds=antic_factor *  dx 
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


