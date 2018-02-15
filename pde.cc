
/********************************************************************
Some explicit two-level finite difference methods for solving
                  _________________________________________
                  |  d_t uk   = -d_x Fk + Sk, i=1,2       |
                  -----------------------------------------

where uk === u_k(x,t) denote the components of the field vector
********************************************************************/


/********************************************************************/
/*  function  McCormack (uses calc_rhs)                             */
/********************************************************************/


void McCormack (double       u1[], double       u2[], 
                const double F1[], const double F2[],
                const double S1[], const double S2[],
                int     choice_BC,
                double  u1left,  double    u2left, 
                double  u1right, double    u2right,
                double  dt, 
                double  D1[], double D2[],
                int     show_McCormack)

  /* One integration step (time interval dt) of the equation
                  _________________________________________
                  |  d_t uk   = -d_x Fk + Sk, i=1,2       |
                  -----------------------------------------

     The McCormack scheme for left predictor differences
     (choice_method==0) is (omitting k)

     ----------------------------------------------------------------------
     u[i]_pred = u[i]                  - dt*(F[i]-F[i-1])/dx + dt*S[i], 
     u[i]_corr = 1/2*u[i] + 1/2 * (u[i]- dt*(F[i+1]-F[i])/dx + dt*S[i])_pred
     ----------------------------------------------------------------------

     It is of second consistency order for smooth data,
     and of order 2/3 for jumps (LeVeque);
     It is quasi-linearly stable if the CFL conditions
     are fullfilled everywhere; no guarantee can be made
     against nonlinear instability.

     Input: - Fields uk, fluxes Fk, soures Sk at time t  (k=1,2)
            - Kind of boundary conditions (choice_BC)
            - Values of th fields at the boundary ukleft, ukright (k=1,2)
              (they are dummies if the BC are not of Dirichlet type).
            - timestep dt, 
            - numerical diffusivities D1,D2 (zero in GKT or macrosc. IDM)
            
     Output: the fields uk at t + dt.

     Some implementation of "calc_rhs"
     is needed to calculate the corrector with the 
     rhs of the predicted values.
     */ 

{       
  const double durch_dx = 1./dx;
  int    i;
  int    nx    = (int)(xmax* durch_dx);
  double u1pred[NXMAX+1],  u2pred[NXMAX+1];
  double F1pred[NXMAX+1],  F2pred[NXMAX+1];
  double S1pred[NXMAX+1],  S2pred[NXMAX+1];
  double diff1 [NXMAX+1],   diff2[NXMAX+1];

  // Upwind step as prediktor
  // No external diffusion for GKT, macr. IDM, GKT variant,3phases

  if ( ((choice_model==0)||(choice_model==1))||(choice_model==3) || (choice_model==7)
       || (choice_model==8) || (choice_model==9)){
    for (i=1; i<=nx-1; i++){
      u1pred[i] = u1[i] + dt * (S1[i] - (F1[i]-F1[i-1]) * durch_dx );
      u2pred[i] = u2[i] + dt * (S2[i] - (F2[i]-F2[i-1]) * durch_dx );
    }
  } 

  else{
    for (i=1; i<=nx-1; i++) {          //numerical diffusion 
      diff1[i] = D1[i] * (u1[i-1]-2.*u1[i]+u1[i+1]) / SQR(dx);
      diff2[i] = D2[i] * (u2[i-1]-2.*u2[i]+u2[i+1]) / SQR(dx);
    }

    for (i=1; i<=nx-1; i++){           // upwind step
      u1pred[i] = u1[i] + dt*(- (F1[i]-F1[i-1]) * durch_dx + S1[i] + diff1[i]);
      u2pred[i] = u2[i] + dt*(- (F2[i]-F2[i-1]) * durch_dx + S2[i] + diff2[i]);
    }
  }
  boundary_vals(u1pred,u2pred,nx,choice_BC, u1left, u2left, u1right, u2right);

  // Testcode

  if (show_McCormack) {
    printf("\nMcCormack: Prediktor\n"); 
    show_advance(u1, u2, u1pred, u2pred, F1, F2, S1, S2);
  }

  // Calculate rhs. for the predicted values with upwind differences 
  // (downwind = false in calc_rhs) 

  calc_rhs (choice_model, false, u1pred, u2pred, F1pred,F2pred, 
                             S1pred, S2pred, show_McCormack);

  // Corrector= 1/2(Predictor) + 1/2(downwind step with predicted values)

  if ( ((choice_model==0)||(choice_model==1))||(choice_model==3) || (choice_model==7) || (choice_model==9) ){
    for (i=1; i<=nx-1; i++)
    {
      u1[i] = 0.5 * (u1pred[i] + u1[i]
        + dt * (S1pred[i] - (F1pred[i+1]-F1pred[i])* durch_dx)  );
      u2[i] = 0.5 * (u2pred[i] + u2[i]
        + dt * (S2pred[i] - (F2pred[i+1]-F2pred[i])* durch_dx)  );
    }
  }

  else{
    for (i=1; i<=nx-1; i++)
    {
      diff1[i] = D1[i] * (u1pred[i-1]-2.*u1pred[i]+u1pred[i+1]) / SQR(dx);
      diff2[i] = D2[i] * (u2pred[i-1]-2.*u2pred[i]+u2pred[i+1]) / SQR(dx);
    }
    for (i=1; i<=nx-1; i++)
    {
      u1[i] = 0.5 * (u1pred[i] + u1[i]
        + dt * (S1pred[i] - (F1pred[i+1]-F1pred[i])* durch_dx + diff1[i]));
      u2[i] = 0.5 * (u2pred[i] + u2[i]
        + dt * (S2pred[i] - (F2pred[i+1]-F2pred[i])* durch_dx + diff2[i]));
    }
  }

  boundary_vals(u1,u2,nx,choice_BC, u1left, u2left, u1right, u2right);


  // Testcode 

  if (show_McCormack) {
    printf("\nMcCormack: Nach Korrektor\n"); 
    show_advance(u1pred,u2pred,u1,u2,
                 F1,F2,S1,S2);
  }
 
}





void Upwind (
                double       u1[], double       u2[],  
                const double F1[], const double F2[],
                const double S1[], const double S2[],
                int choice_BC,
                double     u1left, double     u2left, 
                double    u1right, double    u2right, 
                double dt,
                double D1[], double D2[],
                int    show)

  /*
  The upwind scheme for updating the fields uk is (omitting k)

      ---------------------------------------------------
      (u[i])_new = u[i] - dt (F[i] - F[i-1])/dx + dt*S[i] 
      ---------------------------------------------------

  consistency order: 1 for smooth data, 1/2 for shocks (LeVeque)
  Stability criteria:
    (i)   The eigenvalues lambda of the functional matrix 
          d(F1,F2) / d(u1,u2) must all be positive
    (ii)  Convective CFL conditions |lambda_max| dt/dx < 1
    (iii) Diffusive CFL condition Dmax dt/(2dx^2) < 1
  Terminology and variables explained in the head of "McCormack".
  */

{
  const double durch_dx = 1./dx;
  int    i;
  int    nx    = (int)(xmax* durch_dx);
  double diff1 [NXMAX+1],   diff2[NXMAX+1];
  double u1start[NXMAX+1],  u2start[NXMAX+1];

  if (show) cprhoQ(u1,u2,u1start,u2start);

  // No diffusion for GKT, macr. IDM, GKT variant
  // !! set "choice_model==0" as first condition if flow-violating
  // high-dens corr for GKT is used (implemented in timestep)

  if ( ((choice_model==0)||(choice_model==1))||(choice_model==3) || (choice_model==7)
       || (choice_model==8) || (choice_model==9)){
    for (i=1; i<=nx-1; i++){
      u1[i] = u1[i] + dt * (S1[i] - (F1[i]-F1[i-1]) *durch_dx );
      u2[i] = u2[i] + dt * (S2[i] - (F2[i]-F2[i-1]) *durch_dx );
    }
  }

  // include numerical diffusion for other models
  else{
    for (i=1; i<=nx-1; i++) {
      diff1[i] = D1[i] * (u1[i-1]-2.*u1[i]+u1[i+1]) / SQR(dx);
      diff2[i] = D2[i] * (u2[i-1]-2.*u2[i]+u2[i+1]) / SQR(dx);
    }

    for (i=1; i<=nx-1; i++){
      u1[i] = u1[i] + dt*(- (F1[i]-F1[i-1]) *durch_dx + S1[i] + diff1[i]);
      u2[i] = u2[i] + dt*(- (F2[i]-F2[i-1]) *durch_dx + S2[i] + diff2[i]);
    }
  }

  boundary_vals(u1,u2,nx,choice_BC, u1left, u2left, u1right, u2right);

  if (show) 
  {
    printf("\n Upwind: \n"); 
    show_advance(u1start,u2start,u1,u2,
                 F1,F2,S1,S2);
  }

}


//*******************************************************



void show_advance(const double u1start[], const double u2start[],
                  const double u1step[], const double u2step[], 
                  const double F1[], const double F2[], 
                  const double S1[], const double S2[])


    // uses always left differences for first derivatives 

  {

    int i;

    printf("\nbefore advance by dt: \n");

    printf("ix\t1000rho\t Q\t  v \n");

    for (i=ix_show_min;i<=ix_show_max;i++)
      printf("%i\t %.3f\t %.4f\t  %.4f \n",
              i, 
              1000.*u1start[i], 
              u2start[i], 
              u2start[i]/u1start[i]); 

    printf("\nafter advance by dt: \n");

    printf("\n ix\t   1000*rho\t   Q\t  v \n");

    for (i=ix_show_min;i<=ix_show_max;i++)
      printf("%i\t %.3f\t   %.4f\t  %.4f \n",
              i, 
              1000.*u1step[i], 
              u2step[i], 
              u2step[i]/u1step[i]);

    printf("\n ix\t -1000*d_x F1  -d_x F2    S1   S2   dt*sum2  \n");

    int imindiff = (ix_show_min==0) ? 1 : ix_show_min;
    for (i=imindiff;i<=ix_show_max;i++)
      printf("%i   %.9f  %.5f  %.5f  %.5f    %.5f  \n",
      i, 
      - 1000.*(F1[i]-F1[i-1])/dx, 
      - (F2[i]-F2[i-1])/dx, 
      S1[i],
      S2[i], 
      u2step[i]-u2start[i]); 

  }

