

// *******************************************************************
//  function  input (needs filecheck, error from general.cc) 
// ********************************************************************

 #include <cstring>
 #include <cstdlib>


using namespace std;

int input(char projName[] )
{
  FILE * fp;                    
  char   param_file[MAXSTR];         // Name of the parameter file 
  int i;

  sprintf(param_file,"%s.inp",projName);
  printf("Reading simulation parameters from file %s\n",param_file);

  fp=fopen(param_file,"r");
  filecheck(fp,param_file);

  // GKT model parameters 

  getvar(fp,&v0);
  getvar(fp,&Tr);
  getvar(fp,&tau0);
  getvar(fp,&A0);
  getvar(fp,&rhomax);
  getvar(fp,&antic_factor);
  getvar(fp,&dA);

  // initial and boundary conditions, inflow 

  getvar(fp,&rho_init);
  getvar(fp,&ampl_init);
  getvar(fp,&w_init);
  getvar(fp,&wavel_init);
  getvar(fp,&choice_init);
  getvar(fp,&center_init);
  getvar(fp,&choice_BC);
  getvar(fp,&choice_rmp); 

// Range of integration and numerical parameters 

  getvar(fp,&tmax);
  getvar(fp,&dt);
  getvar(fp,&xmax);
  getvar(fp,&dx);
  getvar(fp,&D1);
  getvar(fp,&D2);

// control of the output (what; in which form) 

  getvar(fp,&dtout); dtout = dtout*(1+1/HUGE_VALUE);
  getvar(fp,&dxout); dxout = dxout*(1+1/HUGE_VALUE);
  getvar(fp,&choice_outp);

  double dtloop_preferred=60;
  dtloop = (tmax/dtloop_preferred<NTLOOPMAX) ? dtloop_preferred : tmax/NTLOOPMAX;

// select variants of the model and initial conditions/BC 

  getvar(fp,&choice_model);  
  getvar(fp,&choice_method);
  getvar(fp,&choice_A);  

// Flags for testing 

  getvar(fp,&test_timestep);  
  getvar(fp,&test_shift);

// Merging Length for older input files (choice_rmp=1)

  dx_rmp=200.0;                       
  if ((choice_rmp==1)&&(!feof(fp)))  getvar(fp,&dx_rmp);


  fclose(fp);


  // **************************************************************
  // crude test of input 
  // (many other not tested combinations can lead to crash!)
  // **************************************************************

  cout << endl;

     // circle -> wavelength must be multiple of circumference xmax
  if (choice_BC==0) wavel_init = wavel_init * (xmax-dx)/xmax;

  if (antic_factor<0) error
     ("Error: anticipation factor must be >=0 ");
  //  if (dt>1) error 
  //     ("Error: discretisation time dt must be <=1(s) for useful results");
  if (dx/dt < v0) error 
     ("Error: Convective Courant criterium: dx/dt > v0 is not satisfied. ");
  if (SQR(dx)/dt < sqrt(SQR(D1)+SQR(D2)) ) error 
    ("Error: Diffusive Courant criterium: dx^2/dt < num diffusivities not satisfied");
  if ( (int)(xmax/dx) > NXMAX) error
     ("Error: more than NXMAX x segments  not possible. Increase const NXMAX");

  if(choice_outp ==1) 
  {
    if (n_xcuts>NXCUTMAX) error
      ("Error: more than NXCUTMAX+1 x cuts not possible");
  }


  // ***************************************************************
  // get coordinates xmin (upstream boundary), tmin (starting time) 
  // (set to zero, if file with extension .xtmin is missing)
  // apply BEFORE get_xdependence and getting ramp data
  // ***************************************************************

  get_xt_start(projName, xmin, tmin);
 

  // ***************************************************************
  //  get additional Parameters in case of non-GKT models            
  // ***************************************************************

  if ((choice_model==3)||(choice_model==7)) idm.getIDMparams(projName);
  if  (choice_model==4) ovm.getOVMparams(projName);  
  if  (choice_model==2) kkl.getKKLparams(projName);  
  if  (choice_model==8) cout <<"io.cc: Mac3phases reads param by itself"<<endl;

  // **************************************************************
  // get boundary conditions (nothing needed if choice_BC==0, 4, or 6)
  // **************************************************************

  // Inflow: rho Dirichlet, Q in equilibrium. Outflow: hom VN or free

  if ( (choice_BC==1) || (choice_BC==2)) 
  {
     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.BCl",projName);
     get_array(in_fname, n_jumps_l, times_l, rhoBC_l); 
  }

  // Inflow and outflow: rho Dirichlet, Q in equilibrium

  else if (choice_BC==3) 
  {
     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.BCl",projName);
     get_array(in_fname, n_jumps_l, times_l, rhoBC_l);
     sprintf(in_fname,"%s.BCr",projName);
     get_array(in_fname, n_jumps_r, times_r, rhoBC_r);
   }

  // Inflow: Dirichlet for rho and Q, Outflow: Homog. Von-Neumann

  else if (choice_BC == 7)
  {
     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.BCl",projName);
     get_array(in_fname, n_jumps_l, times_l, rhoBC_l , QBC_l);
  }


  // Inflow and outflow: Boundary info for both rho and Q needed

  else if ((choice_BC == 5) || (choice_BC == 8) 
        || (choice_BC == 9) || (choice_BC == 10))
  {
     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.BCl",projName);
     get_array(in_fname, n_jumps_l, times_l, rhoBC_l , QBC_l);
     sprintf(in_fname,"%s.BCr",projName);
     get_array(in_fname, n_jumps_r, times_r, rhoBC_r , QBC_r);
  }


  if ( (choice_BC<0) || (choice_BC>10))
      error("input.cc: Error:  choice_BC<0 and >9 not implemented!");


  // **************************************************************
  // Convert rho (1/km) and Q(1/h) to SI units
  // For backward compatibility, convert only input values 
  // satisfying rho>rholimit and Q>Qlimit. 
  // take first value as criterion
  // Then test BC for the right ranges
  // **************************************************************


     const double rholimit=0.5;   // if rho>rholimit: units 1/km assumed
     const double Qlimit  =1.;    // if Q>Qlimit: units 1/h assumed
     bool rho_in_invkm    = false;
     bool Q_in_invh       = false;

     // detemine if boundaries are in SI units

     bool noTrafo=( (choice_model!=5)||(choice_model<50));
     for(int j=0; j<=n_jumps_l; j++){
       if ((rhoBC_l[j]>=rholimit)&&noTrafo) rho_in_invkm = true;
       if (  QBC_l[j]>=Qlimit)   Q_in_invh = true;
     }
     for(int j=0; j<=n_jumps_r; j++){
       if ((rhoBC_l[j]>=rholimit)&&noTrafo) rho_in_invkm = true;
       if (  QBC_r[j]>=Qlimit)   Q_in_invh = true;
     }

     if (rho_in_invkm){
        for (i=0; i<= n_jumps_l; i++) rhoBC_l[i] /=1000.;
        for (i=0; i<= n_jumps_r; i++) rhoBC_r[i] /=1000.;
     }
     if (Q_in_invh){
        for (i=0; i<= n_jumps_l; i++) QBC_l[i] /=3600.;
        for (i=0; i<= n_jumps_r; i++) QBC_r[i] /=3600.;
     }

     // redefine BC condition 4 ("f-f") only for black scholes
     if(choice_model==5){

       // fetch IC (only file IC choice_init-3, sensible here)
       double rhoICl=0;
       double rhoICr=0;

       if(choice_init==3){
         int    n_IC;              
         double pos[NBC], rhoIC[NBC];
         char   in_fname[MAXSTR];
	 sprintf(in_fname,"%s.IC",projName);
	 get_array(in_fname, n_IC, pos, rhoIC);
         if((n_IC-1)>=NBC){
           error("Error: too many IC data in .IC file (more than NBC lines)");
	 }
         rhoICl = intpextp (pos,rhoIC,n_IC+1, 0); //!!! n_jumps+1
         rhoICr = intpextp (pos,rhoIC,n_IC+1, xmax); //!!! n_jumps+1

       }


       double r0=Tr;
       double xBasis_min=tau0;
       double xBasis_max=xBasis_min + xmax;

       int nbc=NBC/10;

       if(choice_BC==6){

       // left boundary, assume linear "put" profile
         
	 n_jumps_l=nbc-1; // since old convention for loop up to "<="
         for (int i=0; i<nbc; i++){
	   double t=tmax*i/(double)nbc;
	   times_l[i]=t;
	   rhoBC_l[i]=rhoICl*exp(-r0*t);
	 }
       }

       if(choice_BC==6){
         // right boundary for linear "call" profile
	 n_jumps_r=nbc-1; // since old convention for loop up to "<="
	 double xs=xBasis_max-rhoICr;
         for (int i=0; i<nbc; i++){
	   double t=tmax*i/(double)nbc;
	   times_r[i]=t;
	   rhoBC_r[i]=xBasis_max - xs*exp(-r0*t);
	   cout <<"i="<<i<<" rhoBC_r[i]="<<rhoBC_r[i]<<endl;
	 }
	 cout <<"xs="<<xs<<" xBasis_max="<<xBasis_max<<endl;
       }

       // redefine effective BC from free to Dirichlet
       if(choice_BC==6){
         choice_BC=3;
       }


       //exit(0);
     }


     // check input for inconsistencies

     for (i=0; i<= n_jumps_l; i++){
        if (rhoBC_l[i]>=rhomax) error(" rho at inflow >= rhomax!");
        if (rhoBC_l[i]< TINY_VALUE) rhoBC_l[i]=TINY_VALUE;
        if (QBC_l[i]>=Qlimit) error(" Q at inflow >= Qlimit!");
        if (QBC_l[i]<v0*TINY_VALUE) QBC_l[i] = v0*TINY_VALUE;
        if (QBC_l[i]/rhoBC_l[i] > 2.*v0) error("v=Q/rho at inflow > 2*v0");
     }
     for (i=0; i<= n_jumps_r; i++){
        if (rhoBC_r[i]>=rhomax) error(" rho at outflow >= rhomax!");
        if (rhoBC_r[i]< TINY_VALUE) rhoBC_r[i]=TINY_VALUE;
        if (QBC_r[i]>=Qlimit) error(" Q at outflow >= Qlimit!");
        if (QBC_r[i]<v0*TINY_VALUE) QBC_r[i] = v0*TINY_VALUE;
     }
     // calculate the downstream boundary condition for the velocity

     for (i=0; i<= n_jumps_r; i++){
        vBC_r[i] = (rhoBC_r[i]>TINY_VALUE) ? QBC_r[i]/rhoBC_r[i] : v0;
        if (vBC_r[i] >= v0 - TINY_VALUE) vBC_r[i] = v0 - TINY_VALUE;
     }

  // **************************************************************
  // if files <projName>.v0, <projName>.Tr exist,
  // get model-parameter gradients of v0 or T representing flow-conserving 
  // bottlenecks. 
  // if files <projName>.Ts and <projName>.T1, <projName>.T2, ... exist,
  // override possible info in <projName>.Tr
  // and get time-dependend variations of Tr (third line)
  // If none of the above applies, use
  // the the constant parameters from <projName>.inp
  // (which are overriddenotherwise).
  // **************************************************************

  char v0str[3]="v0", Trstr[3]="Tr";
  get_xdependence(projName, v0str, v0, v0_loc);
  get_xdependence(projName,  Trstr, Tr, Tr_loc);
  get_T_xtdependence (projName, Tr, Tr_loc,flag_Txt);


  // **************************************************************
  // get dynamic regulation of the desired velocity if file 
  // <projName>.dyn_v0 exists
  // **************************************************************

  get_dyn_v0(projName, x_v0, v0_x, delta_v0_x,thrsh_x);


  // **************************************************************
  // get time dependent truck fraction
  // from file <projName>.fracTruck, and truck model parameters
  // from file <projName>.parTruck, if these files exist.
  // Otherwise assume 0 percent trucks
  // (In any case, the car model parameters are taken from <projName>.inp).
  // **************************************************************

  trucks.getTruckFrac(projName);
  if(trucks.varTruck == true) 
    trucks.getTruckPar(projName);

  // **************************************************************
  // get desired time slices and space cuts for output (if choice_outp==1);
  // **************************************************************

  if (choice_outp==1)
  {    
     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.t",projName);
     ifstream  infile (in_fname, ios::in);
     if( !infile){       
        cout <<"No file "<<in_fname <<" detected -> "
             <<"save only initial conditions in " << projName << ".t0"<<endl;
        n_slices=0;
        time_slices[0]=0.;
     }
     else{
        get_array(in_fname, n_slices, time_slices);
        if (n_slices>=NXCUTMAX) error ("number of t slices >= NXCUTMAX");
     }
     sprintf(in_fname,"%s.x",projName);
     ifstream  infile2 (in_fname, ios::in);
     if( !infile2){       
        cout <<"No file "<<in_fname <<" detected -> "
             <<"save only upstream boundary conditions in " 
             << projName << ".x0"<<endl;
        n_xcuts=0;
        x_xcuts[0]=0.;
     }
     else{
       get_array(in_fname, n_xcuts, x_xcuts);
       if (n_xcuts>=NXCUTMAX) error ("number of x cuts >= NXCUTMAX");
     }
  }


  // **************************************************************
  // get ramp data 
  // choice_rmp==0: no ramps; 1: 1 ramp; 2: several ramps
  // choice_rmp==3: ramps with inflow control (data from <projName>.rmps_c)
  // **************************************************************

  if (choice_rmp==1)  // one ramp at xmax/2
    {
      //  x_rmp:        Position of the ramp
      //  n_jumps_rmp:  Number of data of the ramp-flow time series
      //  times_rmp[i]: the corresponding time values
      //  Q_rmp[i]:     the flow on the ramp for t in veh/h

     x_rmp  = xmax/2.;
     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.rmp",projName);
     get_array(in_fname, n_jumps_rmp, times_rmp, Q_rmp);
     if(n_jumps_rmp>NBC) error ("to many (>NBC) different ramp levels!");
     for(i=0; i<= n_jumps_rmp; i++) Q_rmp[i] /=3600.;
    }


  if ((choice_rmp==2)||(choice_rmp==3))
  {
      // n_rmps:            Number of ramps
      // x_rmps      []:    Positions of the ramps (centered)
      // dx_rmps     []:    Merging lengths of the ramps
      // n_jumps_rmps[]:    Number of data of the ramp-flow time series
      //                    (each element for one ramp)
      // Q_rmps      [][]:  Ramp-flow time series (veh/h) for each ramp
      // times_rmps  [][]:  the corresponding time values for each ramp

     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.rmps",projName);  

     // if controlled ramps, no "uncontrolled ramps" necessary

     ifstream  infile (in_fname, ios::in);
     if((choice_rmp==3)&&(!infile))
       n_rmps=-1;
     else get_array(in_fname, n_rmps, x_rmps, dx_rmps);

     if (n_rmps > NRMPMAX) error ("number of ramps > NRMPMAX");

     // shifting positions to internally calculate with x in [0,xmax]
     int  i_rmps=0;
     for (i_rmps=0; i_rmps<= n_rmps; i_rmps++)  x_rmps[i_rmps]-=xmin;

     // reading ramp flows
     for (i_rmps=0; i_rmps<= n_rmps; i_rmps++)
     {
       sprintf(in_fname,"%s.rmp%i",projName,i_rmps+1);
       get_array(in_fname, n_jumps_rmps[i_rmps], 
                           times_rmps  [i_rmps], 
                           Q_rmps      [i_rmps]);
       if(n_jumps_rmps[i_rmps]>NBC) {
	 cerr <<"to many flow data (>NBC) in time series of ramp "
              <<i_rmps << endl; exit(-1);
       }
     }

     // scaling ramp flows to internally used SI units   

     int  i_time=0;
     for(i_rmps=0; i_rmps<= n_rmps; i_rmps++)
     {
       for(i_time=0; i_time<= n_jumps_rmps[i_rmps]; i_time++) 
          Q_rmps[i_rmps][i_time] /=3600.;
     }
  }


 if (choice_rmp==3)   
    //Several controlled ramps from file <projName>.rmps_c
   {
      // n_rmps_c:            Number of ramps
      // x_rmps_c      []:    Positions of the ramps
      // dx_rmps_c     []:    Merging lengths of the ramps
      // n_jumps_rmps_c[]:    Number of data of the ramp-flow time series
      // Q_rmps_c      [][]:  Ramp-flow time series (veh/h) for each ramp
      // times_rmps_c  [][]:  the corresponding time values for each ramp

     char in_fname[MAXSTR];
     sprintf(in_fname,"%s.rmps_c",projName);
     get_array(in_fname, n_rmps_c, x_rmps_c, dx_rmps_c, s1_rmps_c, max_flow);
     if (n_rmps_c > NRMPMAX) error ("number of ramps > NRMPMAX");
   
     int  i_rmps=0;
     for (i_rmps=0; i_rmps<= n_rmps_c; i_rmps++)
       {
       sprintf(in_fname,"%s.rmp_c%i",projName,i_rmps+1);
       get_array(in_fname, n_jumps_rmps_c[i_rmps], 
                           times_rmps_c  [i_rmps], 
		 Q_rmps_c      [i_rmps]); 

       if(n_jumps_rmps_c[i_rmps]>NBC) {
         cerr <<"to many flow data (>NBC) in time series of controlled ramp "
              <<i_rmps << endl; exit(-1);
         }

       }
    
     int  i_time=0;
     for(i_rmps=0; i_rmps<= n_rmps_c; i_rmps++)
     {
       for(i_time=0; i_time<= n_jumps_rmps_c[i_rmps]; i_time++) 
          Q_rmps_c[i_rmps][i_time] /=3600.;
     }
  }


 
  // *************************************************************
  // Interactive input for debugging if switch test_timestep on
  // **************************************************************

  if(test_timestep)
  {
    cout << "enter it_show_min,it_show_max,ix_show_min,ix_show_max" <<endl
         << "("
         << " it = 0 ... " << (int) (tmax/dt)
         << " ix = 0 ... " << (int) (xmax/dx)
         << " )" << endl;
    cin  >> it_show_min >> it_show_max >> ix_show_min >> ix_show_max;
  };

  if(true){
    for (int i=0; i<=n_jumps_r; i++){
      cout <<"i="<<i<<" rhoBC_r[i]="<<rhoBC_r[i]<<endl;
    }
  }

  return(0);
}

// *************************************************************
// end method input
// *************************************************************







 // **************************************************************
 //       function  get_xt_start 
 // **************************************************************


void   get_xt_start(char projName[], double& xmin, double& tmin)
{
   char in_fname[MAXSTR];

   sprintf(in_fname,"%s.%s",projName, "xtmin");
   ifstream  infile (in_fname, ios::in);

   if( !infile) 
   {
     xmin = 0.;
     tmin = 0.;
     cout << "No file "<<in_fname<<" detected -> xmin=tmin=0" << endl;
   }
   else
   {
     int n_tab;
     double min_tab[2];
     get_array (in_fname, n_tab, min_tab);
     xmin = min_tab[0];
     tmin = min_tab[1];
     if(fabs(tmin)>TINY_VALUE) {
       tmin = 0.; 
       cout << "variable tmin not yet implemented; set to zero" << endl;
     }
     cout << endl <<"NOTICE: read file " << in_fname 
          << " containing xmin = " << xmin 
          << ", tmin = " << tmin << endl << endl;
   }
}



// **************************************************************
//       function  get_xdependence
// **************************************************************

// get x-dependent model parameters v0, Tr (etc)
// from files <projName>.Tr, <projName>.v0 etc,
// if these files exist.
// Otherwise use the constant parameter var_inp

void   get_xdependence(char projName[], char varname[], 
                     double var_inp, double var_field[])
{
   char in_fname[MAXSTR];
   int  ix;

   sprintf(in_fname,"%s.%s",projName, varname);
   ifstream  infile (in_fname, ios::in);

   if( !infile)             // no .Tr file -> use constant var_inp from  .inp
   {
     for (ix=0; ix<=(int)(xmax/dx); ix++) var_field[ix] = var_inp;
     cout << "No file "<<in_fname<<" detected";
     if (varname[0]!='T')
       cout << " -> " << varname << " = const = " << var_inp << endl;
     else cout << endl;
   }
   else
   {
     int n_tab;
     double x_tab[NXCUTMAX+1], var_tab[NXCUTMAX+1];
     get_array (in_fname, n_tab, x_tab, var_tab);
     for (ix=0; ix<=(int)(xmax/dx); ix++) 
       var_field[ix] = intpextp(x_tab,var_tab,n_tab+1,ix*dx+xmin);

     cout << endl <<"NOTICE: read file " << in_fname 
          << " containing an x dependent " 
          << "model parameter:" << endl 
          << " " << varname << " = variable; " 
          << varname << "_loc[0] = " << var_field[0] << " "
          << varname << "_loc[nx] = " << var_field[(int)(xmax/dx)] 
          << endl;
   }
}


 //############################################################
 //       function  get_T_xtdependence
 //############################################################

// get x AND t -dependent model parameter Tr
// from files <projName>.Ts, <projName>.T1 (analog to .rmps, .rmp1..)
// if these files exist.
// Otherwise, use either time-independent parameters if 
// <projName>.Tr exists [via get_xdependence],
// or the constant value from <projName>.inp) 

void get_T_xtdependence (char projName[], double Tr_inp, 
                       double Tr_loc[], bool& flag_Txt)

{
     //n_Ts;                      // how many x sections with different Tr
     //x_Ts[NRMPMAX];             // center of these x sections
     //dx_Ts[NRMPMAX];            // length of these x sections
     //width_Ts[NRMPMAX];         // transition width of slope of Tr
     //n_jumps_Ts[NRMPMAX];       // number of jumps in the sections
     //times_Ts[NRMPMAX][NBC];    // times of jumps in the sections
     //dT_Ts[NRMPMAX][NBC];       // differenc of Tr in .inp

   char in_fname[MAXSTR];
   sprintf(in_fname,"%s.Ts",projName);
   ifstream  infile (in_fname, ios::in);

   if( !infile)            
   {
     flag_Txt=false;
     cout << "No file "<<in_fname<<" detected -> Tr not time dependent"
          << endl;
   }

   else
   {
     cout << endl << endl << "NOTICE: File "<<in_fname
          << " detected: Tr depends on time and location!" << endl;
     flag_Txt=true;
     get_array(in_fname, n_Ts, x_Ts, dx_Ts, width_Ts);
     if (n_Ts > NRMPMAX) error ("number of variable Tr sections > NRMPMAX");

     // shifting positions to internally calculate with x in [0,xmax]
     int  i_Ts=0;
     for (i_Ts=0; i_Ts<= n_Ts; i_Ts++)  x_Ts[i_Ts]-=xmin;

     // reading ramp flows
     for (i_Ts=0; i_Ts<= n_Ts; i_Ts++)
     {
       sprintf(in_fname,"%s.T%i",projName,i_Ts+1);
       get_array(in_fname, n_jumps_Ts[i_Ts], 
                           times_Ts  [i_Ts], 
                           dT_Ts      [i_Ts]);
     }
   
     for(i_Ts=0; i_Ts<= n_Ts; i_Ts++)
     {
       if(n_jumps_Ts[i_Ts]>NBC) 
          error ("to many (>NBC) different T values!");
     }
   }  // end else
} // end get_T_xtdependence


//############################################################
// Get velocity-control data (if <projName>.v0_dyn exists)
//############################################################

void get_dyn_v0(char projName[], double x_v0[], double v0_x[],
		double delta_v0_x[],double thrsh_x[])
{
char in_fname[MAXSTR];
sprintf(in_fname,"%s.v0_dyn",projName);
   ifstream  infile (in_fname, ios::in);

   if( !infile) 
     {
       printf("No dynamic regulation of the desired velocity\n");
       choice_dyn_v0=0;
     }
   else
     {
       get_array (in_fname, n_jmps_v0, x_v0, v0_x, delta_v0_x, thrsh_x);
       choice_dyn_v0=1;
     }
}


//############################################################
//  function getvar
//############################################################

// interprets lines beginning with COMMENTCHAR as comments;
// reads the first word of other lines;
// inteprets the remaining contents also of these lines as comments

void getvar(FILE *fp, double *pdouble)
{
  char   comment[MAXSTR];               // Dummy string for comments 
  if (fgets(comment,MAXSTR,fp))
  {
    while ((comment[0]==COMMENTCHAR) && fgets(comment,MAXSTR,fp)){;}
    sscanf(comment,"%lf",pdouble);
    printf("%.3f\t %s",*pdouble,comment);
  }
}


void getvar(FILE *fp, int *pint)
{
  char   comment[MAXSTR];    
  if (fgets(comment,MAXSTR,fp))
  {
    while ((comment[0]==COMMENTCHAR) &&fgets(comment,MAXSTR,fp)){;}
    sscanf(comment,"%i",pint);
    printf("%i\t %s",*pint,comment);
  }
}




// ********************************************************************
//  function get_array                      
// ********************************************************************

// Opens file "fname" and defines with its data "array"
//   or two arrays "array1","array2", and the length n_array+1
//   = number of lines of the file
//
//   The format of the file must be 
//   value0
//   value1
//   ...
//   or
//   value01 value02
//   value11 value12
//   ...
//
//   Comment lines begining with COMMENTCHAR='%' are allowed
//


int is_not_white(char line[], int nchar)
{
  int not_white=false;
  for (int i=0; i<=nchar-1; i++)  // last always '\n'
  {
    //  cout << " line["<<i<<"]="<<line[i] << "nonwhite = "
    //      <<( ( (line[i]!=' ') && (line[i] !='\n') )? 1 : 0) << endl;
    if( (line[i]!=' ') && (line[i] !='\n') ) not_white=true;
  }
  return(not_white);
}


void  get_array (char* fname, int & n_array, double array[])
{
  cout << endl << "in get_array: Read file "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[MAXSTR];
  while (infile.getline(line,MAXSTR))  // One input line
  {
  if(is_not_white(line,strlen(line)))
    if( (line[0] != COMMENTCHAR) && (strlen(line)>0) )
    {
      std::istringstream linestream (line);
      linestream >> array[i];
   //    cout << "is_not_white = " << is_not_white(line,strlen(line)) << " ";
      cout << "array[" << i <<"] = " << array[i] << endl;
      i++;
    }
  }
  infile.close();
  n_array = i-1;
  // cout << "n_array = " << n_array << endl;
}

void  get_array (char* fname, int & n_array, 
      double array1[], double array2[])
{
  cout << endl << "in get_array: Read file "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[MAXSTR];
  while (infile.getline(line,MAXSTR))  // One input line
  {
  if(is_not_white(line,strlen(line)))
    if( (line[0] != COMMENTCHAR) && (strlen(line)>0) )
    {
      std::istringstream linestream (line);
      linestream >> array1[i] >> array2[i];
      cout << "array1[" << i <<"] = " << array1[i] 
           << " array2[" << i <<"] = " << array2[i] << endl;
      i++;
    }
  }
  infile.close();
  n_array = i-1;
  // cout << "n_array = " << n_array << endl;
}


void  get_array (char* fname, int & n_array, 
      double array1[], double array2[], double array3[])
{
  cout << endl << "in get_array: Read file "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[MAXSTR];
  while (infile.getline(line,MAXSTR))  // One input line
  {
  if(is_not_white(line,strlen(line)))
    if( (line[0] != COMMENTCHAR) && (strlen(line)>0) )
    {
       std::istringstream linestream (line);
       linestream >> array1[i] >> array2[i] >> array3[i];
       cout << "array1[" << i <<"] = " << array1[i] 
            << " array2[" << i <<"] = " << array2[i] 
	    << " array3[" << i <<"] = " << array3[i]<< endl;
      i++;
    }
  }
  infile.close();
  n_array = i-1;
  // cout << "n_array = " << n_array << endl;
}

void  get_array (char* fname, int & n_array, 
      double array1[], double array2[], double array3[], double array4[])
{
  cout << endl << "in get_array: Read file "<<fname <<endl;
  ifstream  infile (fname, ios::in);
  if( !infile)
  {
    cerr << "Error opening file " << fname << " for reading" << endl;
    exit(-1);
  }
  
  int i=0;
  char line[MAXSTR];
  while (infile.getline(line,MAXSTR))  // One input line
  {
  if(is_not_white(line,strlen(line)))
    if( (line[0] != COMMENTCHAR) && (strlen(line)>0) )
    {
       std::istringstream linestream (line);
       linestream >> array1[i] >> array2[i] >> array3[i] >> array4[i];
       cout << "array1[" << i <<"] = " << array1[i] 
            << " array2[" << i <<"] = " << array2[i]
	    << " array3[" << i <<"] = " << array3[i]
	    << " array4[" << i <<"] = " << array4[i]<< endl;
      i++;
    }
  }
  infile.close();
  n_array = i-1;
  //  cout << "n_array = " << n_array << endl;
}


// ##################################################################
//  function write_cross_section                                    
// ##################################################################

// Writes cross-sections of the data at t=const (t_const=true)
//   or x=const (otherwise) to the file 
//   "<projName>.t<t_sec>" or "<projName>.x<x_meter>"
//   If x=const, the field length is tmax/dtloop
//   If t=const, the field length is xmax/dx, and x=x_internal+xmin.
//
//   output  format  (not in SI units!) :
//   x(km), t(min), rho(1/km), v(km/h), Q(veh/h), a(m/s^2)
   

void  write_cross_section (int t_const, double t_sec, double x_meter, 
      const double rho[], const double Q[], const double a[])

{

  bool noTrafo=( (choice_model==5)||(choice_model>50));

  //!! quick hack for Black Scholes (choice_model==5)
  // with other time scales

  dtloop=0.01*tmax;
  int     n_field = (t_const==true) ? (int)(xmax/dx) : (int)(tmax/dtloop);
  //cout << "write_cross_section: t_const="<<t_const<<" xmax="<<xmax<<" dx="<<dx<<" tmax="<<tmax<<" dtloop="<<dtloop<<" n_field = " << n_field << endl;

  char    out_name[MAXSTR];
  double pref_x=(noTrafo)   ? 1 : 0.001;  // 5=Black Scholes or FPE
  double pref_t=(noTrafo)   ? 1 : 1./60.; // 5=Black Scholes or FPE
  double pref_rho=(noTrafo) ? 1 : 1000.;  // 5=Black Scholes or FPE

  if (t_const==true) 
    sprintf(out_name,"%s%s%i",projName, ".t", (int) t_sec);
  else 
    sprintf(out_name,"%s%s%i",projName, ".x", (int) x_meter);


  cout << "in write_cross_section: "
       << "Write to file "<<out_name<<" noTrafo="<<noTrafo<<endl;
  ofstream  outfile (out_name, ios::out);
  if( !outfile)
  {
    cerr << "Error opening file" << out_name << "for writing" << endl;
    exit(-1);
  }


 // title line

  outfile << "# x(km)"<<'\t'<<"t(min)\t\t rho(1/km) \t"
          << "v(km/h)  \t Q(1/h) \t a(m/s^2)" 
          << endl;                       

 // actual writing

  outfile.precision(4);


  for (int i=0; i<n_field; i++)
  {
    double v_kmh        = (rho[i]==0) ? 0 : 3.6 * Q[i]/rho[i];
    double x_km         = (t_const)   ? pref_x*(i*dx) : pref_x*x_meter;
    //                  x_km += pref_x*xmin;
    double t_min        = (t_const)   ? pref_t* t_sec : pref_t*i*dtloop;
    outfile <<  x_km +pref_x*xmin  << "\t" 
            <<  t_min            << "\t\t" 
            <<  pref_rho* rho[i]    << "\t\t" 
            <<  v_kmh            << " \t\t"
            <<  3600.*Q[i]       << " \t\t"
            <<  a[i]       << endl;
  }

  outfile.close();

}



// ##################################################################
//  function  write_results (needs filecheck)                       
// ##################################################################

void write_results(char projName[], int it,
     const double rho[], const double Q[], const double a[])

// ##########################################################
//   Writes (it=0) or appends (it!=0) the results of the last iteration
//     to the file "<projName>.dat";
//     if choice_BC==10 ( comparison with measured data), the additional
//     files "<projName>.left", "<projName>.right", and "<projName>.test"
//     are written
//     output  format  (not in SI units!) :
//     x(km), t(min), rho(1/km), v(km/h), Q(veh/h), a(m/s^2)
// ##########################################################

{
  bool noTrafo=( (choice_model==5));

  FILE  *fp;                  // data(x,t) 
  char   param_file[MAXSTR];
  int    i,ix,ndxout,nxout;
  double pref_x=(noTrafo)   ? 1 : 0.001;  // 5=Black Scholes or FPE
  double pref_t=(noTrafo)   ? 1 : 1./60.; // 5=Black Scholes or FPE
  double pref_rho=(noTrafo) ? 1 : 1000.;  // 5=Black Scholes or FPE

  sprintf(param_file,"%s.dat",projName);


  //  printf("in write_results: data appended to %s\n",param_file); 
  fp     = (it==0) ? fopen(param_file,"w") : fopen(param_file,"a");
  filecheck (fp,param_file);
  ndxout = ((dxout/dx)>=1) ? (int) (dxout/dx) : 1;
  nxout  = (int) (xmax/(ndxout*dx));

  if(it==0) 
    fprintf (fp, 
      "# x(km) \t t(min) \t rho(1/km) \t v(km/h) \t Q(1/h) \t a(m^2/s^2) \n");

  for (i=0;i<=nxout;i++)
  {
    ix            = i*ndxout;
    double v_kmh  = (rho[ix]==0) ? 0 : 3.6 * Q[ix]/rho[ix];

    fprintf(fp, "%.6f\t  %.6f\t  %f\t %f\t %f\t %f\n",
          pref_x*(ix*dx +xmin), 
          pref_t * it*dt, 
          pref_rho* rho[ix], 
          v_kmh,
          3600.*Q[ix],
          a[ix]
          );
  }
  fprintf(fp, "\n");  // newline means new t step for gnuplot 
  fclose(fp);
}

