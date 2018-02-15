


//########################################################################
// trucks.cc: deals with time-dependent fraction of trucks
// leading to time-dependent macroscopic model parameters
//#########################################################################



Trucks::Trucks()
{
  v0TruckCorr=TrTruckCorr=rhomaxTruckCorr=1.;
}


void Trucks::getTruckFrac(char namepar[])
{
  char   fname[MAXSTR];         // Name of the parameter file 
  double dummy2[NBC+1],dummy3[NBC+1];

  sprintf(fname,"%s.fracTruck",namepar);
  ifstream  infile (fname, ios::in);
  if( !infile)                // no .fracTruck file -> use time-indep. v0 etc
  {
    varTruck = false;
    cout <<"No file "<<fname
         <<" detected -> Truck fraction not time dependent (zero)" << endl;
  }
  else
  {
    varTruck = true;
    get_array(fname,nArr,timesArr,dummy2,dummy3,fracTruckArr);
    cout <<"File "<<fname
         <<" detected -> read time-dependent truck fractions from col 4" 
         << endl;
  }
}


void Trucks::getTruckPar(char namepar[])
{
  FILE * fp;                    
  char   param_file[MAXSTR];         // Name of the parameter file 
  sprintf(param_file,"%s.parTruck",namepar);
  printf("Reading fractions of trucks from file %s\n",param_file);
  fp=fopen(param_file,"r");
  filecheck(fp,param_file);

  getvar(fp,&v0Truck);
  getvar(fp,&TrTruck);
  getvar(fp,&rhomaxTruck);
  fclose(fp);
}

// calculate correction factors with respect to pure care traffic

void Trucks::calc_corrections(double t)
{
  double fracTruck = intpextp(timesArr, fracTruckArr, nArr+1, t);
  v0TruckCorr = (fracTruck*v0Truck + (1-fracTruck)*v0) / v0;
  TrTruckCorr = (fracTruck*TrTruck + (1-fracTruck)*Tr) / Tr;
  rhomaxTruckCorr = (fracTruck*rhomaxTruck + (1-fracTruck)*rhomax) / rhomax;
}
