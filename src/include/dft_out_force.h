/*
  Prototypes for dft_out_force.c
*/
void calc_force(FILE *fp, double *x,double fac_area);
void sum_rho_wall(double *x, double **Sum_rho);
void force_elec(double *x, double **Sum_dphi_dx);
void find_offset(int el_type,int jdim,int *offset);
double calc_deriv(int idim,int inode0,int flag,int *blocked, double *x, int ilist);
void find_pot_derivs(double *x, double *psi_deriv);
double sum_rho_midplane(double *x);
void integrate_rho_vdash(double *x,double **rho_vdash);
double calc_local_pressure(double *rho);
