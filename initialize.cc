
/********************************************************************/
/*  function initialize                                             */
/********************************************************************/

void initialize(double rho[], double Q[])

  /* Set initial conditions for rho and Q; 
     kind of initialization chosen by "choice_init:
     choice_init=0:  rho = rhoinit + periodic * Gauss envelope;
                     Q   = Qe(rho)
     choice_init=1:  rho = rhoinit + ampl_init if x in [x1,x2],
                           rhoinit otherwise
                           where x1/x2 = xmax*center_init +/- w_init/2
                     Q   = Qe(rho)
     choice_init=2:  rho = rhoinit + asymmetric dipole-like peturbation 
                     (max amplitude ampl_init) as used by Kerner et. al.
                     Q   = Qe(rhoinit)
     choice_init=3:  rho = data from file <projName>.IC
                     Q   = Qe(rho)
     choice_init=4:  rho = data from file <projName>.IC
                     Q   = data from file <projName>.IC

   */
{

  int ix, nx = (int)(xmax/dx);
  
  if(choice_init==0)
  {
    double x_center = xmax*center_init;
    double invsqrt2=1./sqrt(2.);
    for(ix=0; ix<=nx; ix++)
    {
      rho[ix]   = rho_init + exp(-SQR((ix*dx - x_center)*invsqrt2/w_init))
                      * ampl_init * cos(2.*PI*(ix*dx - x_center)/wavel_init);
      Q[ix] = rho[ix] * intp(veqtab, NRHO+1, rho[ix], 0, rhomax); 
    }
  }

  else if(choice_init==1)
  {
    int ix_start = (int) ((xmax*center_init - w_init/2.)/dx -1); 
    int ix_end   = (int) ((xmax*center_init + w_init/2.)/dx);
    if (ix_start<0) ix_start=0;
    if (ix_end>nx) ix_end=nx;

    for(ix=0;        ix< ix_start;ix++) rho[ix] = rho_init; 
    for(ix=ix_start; ix<=ix_end;  ix++) rho[ix] = rho_init + ampl_init; 
    for(ix=ix_end+1;   ix<= nx;     ix++) rho[ix] = rho_init; 
    for(ix=0; ix<=nx;  ix++) 
         Q[ix] = rho[ix] * intp(veqtab, NRHO+1, rho[ix], 0, rhomax); 
  }


  else if(choice_init==2) // perturbation a la Kerner
  {
    double wplus  = 32200./160.;  // L/160 in Eq( 10) with L=32.2 km
    double wminus = 32200./40.;   // L/40  in Eq( 10) with L=32.2 km
    double dw     = 32200./32.;   // 5L/16-11L/32 
                                  // = distance between pos/neg hump
    for(ix=0; ix<=nx;  ix++) 
         rho[ix] = rho_init 
                 + ampl_init / SQR( cosh( (ix*dx-xmax*center_init)/wplus))
                 - wplus/wminus*ampl_init 
                   / SQR( cosh( (ix*dx-xmax*center_init-dw)/wminus));

    for(ix=0; ix<=nx;  ix++) 
      // use this if flow and velocity should be in local equilibrium
     //  Q[ix] = rho[ix] * intp(veqtab, NRHO, rho[ix], 0, rhomax); 
      // use this if velocity should be const -> flow Q=rho*v
         Q[ix] = rho[ix] * intp(veqtab, NRHO+1, rho_init, 0, rhomax); 
  }

  else if(choice_init==3) // initial rho from file; Q in equilibrium 
  {
    int    n_jumps;               // How many data in .IC file
    double pos[NBC], rhoIC[NBC];  // Positions and corresp. initial densities
    char   in_fname[MAXSTR];

    const double rholimit=0.5;   // if max(rho)>rholimit: units 1/km assumed
    bool rho_in_invkm    = false;

    sprintf(in_fname,"%s.IC",projName);
    get_array(in_fname, n_jumps, pos, rhoIC);
    if((n_jumps-1)>=NBC) 
       error("Error: too many IC data in .IC file (more than NBC lines)");

    // convert to internal SI units

    for(int j=0; j<=n_jumps; j++){
      if ((rhoIC[j]>rholimit)&&(choice_model !=5)&&(choice_model <50)){
        rho_in_invkm = true;
      }
    }
    if (rho_in_invkm){
      cout << "Assuming rho given in units of veh/km in file "
           <<in_fname<<endl;
      for(int j=0; j<=n_jumps; j++) rhoIC[j]/=1000.;
    }

    for(ix=0; ix<=nx;  ix++) {
      rho[ix] = intpextp (pos,rhoIC,n_jumps+1, ix*dx); //!!! n_jumps+1
      Q[ix] = rho[ix] * intp(veqtab, NRHO+1, rho[ix], 0, rhomax); 
      if(false){
	cout <<"initialize: n_jumps="<<n_jumps<<" rho["<<ix<<"]="<<rho[ix]<<endl;
      }
    }
  }

  else if(choice_init==4)  // initial rho and Q from file
  {
    int    n_jumps;              
    double pos[NBC];             
    double rhoIC[NBC], QIC[NBC]; 
    char      in_fname[MAXSTR];

    const double rholimit=0.5;   // if max(rho)>rholimit: units 1/km assumed
    const double Qlimit  =1.;    // if max(Q)>Qlimit: units 1/h assumed
    bool rho_in_invkm    = false;
    bool Q_in_invh       = false;

    sprintf(in_fname,"%s.IC",projName);
    get_array(in_fname, n_jumps, pos, rhoIC, QIC);
    if((n_jumps-1)>=NBC)
       error("Error: too many IC data in .IC file (more than NBC lines)");

    // convert to internal SI units

    for(int j=0; j<=n_jumps; j++){
      if (rhoIC[j]>rholimit) rho_in_invkm = true;
      if (  QIC[j]>Qlimit)     Q_in_invh  = true;
    }
    if (rho_in_invkm){
      cout << "Assuming rho given in units of veh/km in file "
           <<in_fname<<endl;
      for(int j=0; j<=n_jumps; j++) rhoIC[j]/=1000.;
    }
    if (Q_in_invh){
      cout << "Assuming Q given in units of veh/h in file "
           <<in_fname<<endl;
      for(int j=0; j<=n_jumps; j++) QIC[j]/=3600.;
    }


    for(ix=0; ix<=nx;  ix++) {
      rho[ix] = intpextp (pos,rhoIC,n_jumps+1, ix*dx);
      Q[ix]   = intpextp (pos, QIC, n_jumps+1, ix*dx);
    }
  }

  else error("Error: Choice_init<0 or >4 not implemented!");


  // density must be finite!

  //  for(ix=0; ix<=nx;  ix++) if( rho[ix]<1.e-6) rho[ix]=1.e-6;

  // test output
 
  if(false)
  { cout << "initial conditions:" << endl;
    for(ix=0; ix<nx;  ix++)
       cout <<"rho["<<ix<<"] = " << rho[ix] 
            <<"Q["  <<ix<<"] = " << Q[ix] 
            << endl;
    //exit(0);
  }
}


