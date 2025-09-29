
// ************************************************************
void timestep ( int it, double rho[], double Q[], double a[])
// **************************************************************

  /*
      one complete iteration for fixed timestep dt; 

     On input, the arrays rho, Q are the starting values,
     and on output, the updated results. 
     In addition, the acceleration a is outputted for saving to a file
     in the main program.

     timestep calls 
       "advance_dt" for the actual time step, and 
       "calc_rhs" to calculate the fluxes and sources.

     "advance_dt", in turn, choses an integration method depending on the
     global variable "choice_method":
     choice_method = 0: "McCormack" with left-differentiated predictor
     choice_method = 1: "Upwind" (upwind differencing scheme)

     "calc_rhs" defines the model selected by the
     global variable "choice_model"

     in "timestep", "calc_rhs" is called with downwind differences
     (only relevant if fluxes have derivatives as in the KK model;
     together with the upwind differencing in "Upwind" or in the predictor of
     "McCormack" diffusions are realied by usual central differences)
   */
// **************************************************************

{

  int    nx = (int) (xmax/dx); 
  double rholeft, Qleft, rhoright, Qright;   // Boundary values

  double F1[NXMAX+1],F2[NXMAX+1];             // fluxes
  double S1[NXMAX+1],S2[NXMAX+1];            // sources
  double D1field[NXMAX+1], D2field[NXMAX+1]; // additional high-density diffus


  // determine, if step should be shown for debugging (sic!)

  //show=(it<20);
  show=false;
  if(show){
    cout <<endl<<endl
	 << "in timestep debugging: it= " << it << endl;
  }

  // determine truck corrections if truck percentage is time dependent

  if(trucks.varTruck) trucks.calc_corrections(it*dt);



  // determine x and t dependent T (=Tr), if flag_Txt == true
  // Tr_loc[] must be recalculated in every timestep in this case.

  if (flag_Txt)
  {

    // (1) interpolation in time

    double T_xvals[4*(NRMPMAX+1)];
    double T_field[4*(NRMPMAX+1)];

    for (int i_T=0; i_T<=n_Ts; i_T++)  // n_Ts, x_Ts etc  global var
    {

      // Difference Tr(t) - Tr for each section
      double dT_t =    
        intpextp (times_Ts[i_T], dT_Ts[i_T], n_jumps_Ts[i_T]+1, it*dt);

      // Make arrays T_xvals, T_field of x dependend Tr
      T_xvals[4*i_T]   = x_Ts[i_T] - dx_Ts[i_T]/2 - width_Ts[i_T]/2;
      T_xvals[4*i_T+1] = x_Ts[i_T] - dx_Ts[i_T]/2 + width_Ts[i_T]/2;
      T_xvals[4*i_T+2] = x_Ts[i_T] + dx_Ts[i_T]/2 - width_Ts[i_T]/2;
      T_xvals[4*i_T+3] = x_Ts[i_T] + dx_Ts[i_T]/2 + width_Ts[i_T]/2;
      T_field[4*i_T]   = Tr;
      T_field[4*i_T+1] = Tr + dT_t;
      T_field[4*i_T+2] = Tr + dT_t;
      T_field[4*i_T+3] = Tr;
    }

    // (2) interpolation in x -> Tr_loc

    for (int ix=0; ix<=nx; ix++){
      Tr_loc[ix] = intpextp(T_xvals, T_field, 4*(n_Ts+1), ix*dx);
    }

    // test output

    if (false) {
      for (int ix=0; ix<=nx; ix++) 
         cout << "Tr_loc["<<ix<<"]="<<Tr_loc[ix]<<endl; 
      for (int ix=0; ix<=4*(n_Ts+1)-1; ix++)
	cout << "T_xvals["<<ix<<"]="<<T_xvals[ix]
             <<" T_field["<<ix<<"]="<<T_field[ix]<<endl;
    }
  } // end if (flag_Txt)



  //    on-and off-ramps:


  if (choice_rmp==0);          // no on-or off-ramp

  else if (choice_rmp==1)      // on ramp
  {
    flow_rmp = intpextp (times_rmp, Q_rmp, n_jumps_rmp+1, it*dt);
  }

  else if((choice_rmp==2)||(choice_rmp==3)||(choice_rmp==4))  // several ramps
  {
    int i_rmp;
    for (i_rmp=0; i_rmp<=n_rmps; i_rmp++)
    {
      flow_rmps[i_rmp] = intpextp (times_rmps[i_rmp], Q_rmps[i_rmp], 
        n_jumps_rmps[i_rmp]+1, it*dt);
    }

    if (choice_rmp==3) {             // ramps with flow control
      double dtout_ctl=3;            // logging every dtout_ctl seconds 
      int ndt_ctl = (int)(dtout_ctl/dt);

      for (i_rmp=0; i_rmp<=n_rmps_c; i_rmp++){
         FILE *stfp;
         char stfile[MAXSTR];
         sprintf(stfile,"%s.store%i",projName, i_rmp+1);
         if (it==1) {
           stfp = fopen(stfile,"w");
           fprintf(stfp,
	     "# time(s)\t n_queued\tramp inflow\tcntrl inflow\t");
           fprintf(stfp, "travtime\tnveh_main\n");
           fclose(stfp);
	 }

         if (it%ndt_ctl==0) {
           stfp = fopen(stfile,"a");

           fprintf(stfp,"%f\t",(it*dt));

           double flow_emp = intpextp(times_rmps_c[i_rmp], Q_rmps_c[i_rmp],
		n_jumps_rmps_c[i_rmp]+1, it*dt);

           double flow_sond = intp(Q, int(xmax/dx),s1_rmps_c[i_rmp] ,0,xmax);
           double dens_sond = intp(rho, int(xmax/dx),s1_rmps_c[i_rmp] ,0,xmax);
           flow_rmps_c[i_rmp] = 
             cont_rmpflow(flow_emp,flow_sond,max_flow[i_rmp],
						dens_sond,&rmp_stor[i_rmp]);

	   // instantaneous travelling time integral(dx)(1/V)

           double travtime=0;
           for (int ix=0; ix<=int(xmax/dx); ix++){
             travtime += (Q[ix]>0) ? dx*rho[ix]/Q[ix] : 0;
	   }

	   // total number of cars on main road; difference to number on
	   // free road indicates "number of effectively waiting veh"

           double nveh_main=0;
           for (int ix=0; ix<=int(xmax/dx); ix++){
             nveh_main += dx*rho[ix];
	   }


	   fprintf(stfp,"%f\t%f\t%f\t%f\t%f",rmp_stor[i_rmp],
		      flow_emp,flow_rmps_c[i_rmp], travtime, nveh_main);
           fprintf(stfp,"\n");
           fclose(stfp);
	 }
       }

    }
  }


   // Determine boundary values from tables.


  if((( choice_BC==0) || (choice_BC==4)) || (choice_BC==6))
    rholeft = rhoright = Qleft   = Qright = 0;  // values not needed

  else if (( choice_BC==1) ||  ( choice_BC==2) ||  ( choice_BC==3))
  {
    // Dirichlet BC for density
    rholeft  = intpextp (times_l,rhoBC_l,n_jumps_l+1,it*dt); 
    rhoright = intpextp (times_r,rhoBC_r,n_jumps_r+1,it*dt);
    // equilibr. flow
    Qleft    = rholeft  * intp(veqtab, NRHO+1,rholeft,0,rhomax);
    Qright   = rhoright * intp(veqtab, NRHO+1,rhoright,0,rhomax);
  }

  // choice_BC  upstream BC                  downstream BC
  // -------------------------------------------------------------
 // 8          Dir_rho,  Dir_Q              Dir_rho, equil_Q
  // 9          equilFree_rho(Q), Dir_Q      equilCong_rho(Q), Dir_Q
  // 10         Dir_rho,  Dir_Q              equil_rho(v), equil_Q(v)

  else if( (choice_BC==5) ||  ( choice_BC==7) || (choice_BC==8)) 
  {
    // Dirichlet BC for DENSITY
    rholeft  = intpextp (times_l,rhoBC_l,n_jumps_l+1,it*dt);
    rhoright = intpextp (times_r,rhoBC_r,n_jumps_r+1,it*dt);
    // upstram bundary: nonequilibr. flow
    Qleft    = intpextp (times_l,QBC_l,n_jumps_l+1,it*dt);
    // downstream bundary: equilibr. flow
    Qright   = rhoright * intp(veqtab, NRHO+1,rhoright,0,rhomax);
  }  


  else if( choice_BC==9)
  {
    // Dirichlet BC for FLOW
    Qleft    = intpextp (times_l,QBC_l,n_jumps_l+1,it*dt);
    Qright   = intpextp (times_r,QBC_r,n_jumps_r+1,it*dt);
    // rho upstream: free branch of inverse of Qe(rho) relation
    rholeft  = intp(rho_freetab,NRHO+1,Qleft, 0, QFACTAB*Qmax);
    // rho downstream: congested branch
    rhoright = intp(rho_congtab,NRHO+1,Qright, 0, QFACTAB*Qmax);
  }  

  else if( choice_BC==10)
  {
    // upstream boundary as choice_BC=8
    rholeft  = intpextp (times_l,rhoBC_l,n_jumps_l+1,it*dt);
    Qleft    = intpextp (times_l,QBC_l,n_jumps_l+1,it*dt);

    // Dirichlet downstream BC for VELOCITY
    double vright  = intpextp (times_r,vBC_r,n_jumps_r+1,it*dt);
    rhoright       = intp(rho_vtab,NRHO+1, vright, 0, v0);
    Qright         = rhoright * intp(veqtab, NRHO+1,rhoright,0,rhomax);
  }  

  else error ("timestep: choice_BC<0 or >10 not implemented");


  // Calculate dynamic (density dependent) velocity control

  if (choice_dyn_v0==1) dyn_v0(rho);



  // Include optional numerical diffusion (enable by "if(true)")

  if (true){
    if ( !(((choice_model==0)||(choice_model==1))||(choice_model==3) 
	   || (choice_model==7) || (choice_model==8) || (choice_model==9)
           ) ){
      // no diffus. for choice_model==8!!
      for(int i=0; i<=nx; i++)
      {
        D1field[i] = D1; 
        D2field[i] = D2; 
      }
    }
  }

  // determine flow-violating but more stable  high density correction:
  // if the more unstable but flow-conserving correction in calc_rhs is used:
  // set "if (false)" here and search for "D_corr" in calc_rhs
  // if the method here is used: 
  // change "upwind" and" McCormack" in pde.cc (look for "diffusion")

  if (false){ //!!!
    if (((choice_model==0)||(choice_model==1))||(choice_model==3)||(choice_model==6)
	||(choice_model==42) || (choice_model==7)
	|| (choice_model==8)|| (choice_model==9)){
      double Dplusmax=200; // high-density correction (or in calc_rhs) 
      for(int i=0; i<=nx; i++)
      {
        double high_dens = intp(high_denstab, NRHO+1, rho[i],0,rhomax);
        D1field[i] += Dplusmax * high_dens; 
        D2field[i] += Dplusmax * high_dens; 
      }
    }
  }


  // #######################################################
  // Actual time step
  // #######################################################


  // Calculate rhs (with downwind differences= true)
  // as preparation for the iteration 

  calc_rhs (choice_model, true, rho,Q,F1,F2,S1,S2,show);
 
  advance_dt(it,choice_method,rho,Q,F1,F2,S1,S2, 
            choice_BC,rholeft,Qleft,rhoright,Qright,
            dt,D1field,D2field,show);

  for(int i=0; i<=nx; i++)
    a[i] = S2[i]/rho[i];     // rho=0 catched by calc_rhs


} //end timestep



void cprhoQ (const double rho[], const double Q[],
  double rhocp[], double Qcp[])

  /* Copies  rho,Q to rhocp,Qcpcp */

{
  int i;
  for (i=0; i<=(int)(xmax/dx); i++)
  {
    rhocp[i] = rho[i];
    Qcp[i]   = Q[i];
  }
}



void advance_dt (int it, int choice_method,
		double       u1[], double       u2[],
                const double F1[], const double F2[],
                const double S1[], const double S2[],
                int choice_BC,
                double    u1left,  double    u2left, 
                double    u1right, double    u2right,
                double dt, 
                double D1[], double D2[], 
		bool show)


  { 
    if(choice_model==10){ // overrules .inp file for SGM
      Godunov (u1, u2, F1, F2, S1, S2, 
	       choice_BC, u1left, u2left, u1right, u2right, it, dt, show);
    }
    else if( choice_method == 0) 
    McCormack ( u1, u2, F1, F2, S1, S2, 
         choice_BC, u1left, u2left, u1right, u2right, dt, D1, D2, show);
    else if(choice_method == 1) 
    Upwind (    u1, u2, F1, F2, S1, S2, 
         choice_BC, u1left, u2left, u1right, u2right, dt, D1, D2, show);
  }
