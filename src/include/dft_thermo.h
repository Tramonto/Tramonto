/*
  Prototypes for dft_thermo.c
*/
void  thermodynamics( char *output_file1, int print_flag);
double calc_hs_properties(double *betamu_hs,double *rho);
double calc_att_properties(double *betamu_att, double *rho);
void calc_charge_correlations_b();
double int_stencil_bulk(int sten_type,int icomp,int jcomp);
double uLJatt_n(double r,int i, int j);
double uLJatt_n_int(double r,int i, int j);
double uLJatt_n_noshift(double r,int i, int j);
double deltaC_MSA(double r,int i, int j);
double deltaC_MSA_int(double r,int i, int j);
void pot_parameters(char *output_file1);
double coexistence();
double dp_drho_hs(double *rho);
double dp_drho_att(double *rho);
double dmu_drho_hs(double *rho);
double dmu_drho_att(double *rho);
void print_thermo(char *output_file1, double betap_hs, 
		  double *betamu_hs);
