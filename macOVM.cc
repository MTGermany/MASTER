OVMparams::OVMparams()
{
  V0 = 7.91;      // 7.91
  // V1 = ;
  c0 = 0.16;      // 0.16
  c1 = 1.57;      // 1.57
}


void OVMparams::getOVMparams(char namepar[])
{
  FILE * fp;                    
  char   param_file[MAXSTR];   
  sprintf(param_file,"%s.OVM",namepar);
  cout << "\nOVM model chosen:" << endl
       << "Reading params of OVM model from file "<< param_file
       << endl;
  fp = fopen(param_file,"r");
  filecheck(fp, param_file);

  getvar(fp,&V0);
 //getvar(fp,&V1);
  getvar(fp,&c0);
  getvar(fp,&c1);

  fclose(fp);
}
