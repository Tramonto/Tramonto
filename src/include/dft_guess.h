/*
  Prototypes from dft_guess.c
*/

void set_initial_guess (int iguess,double *x,int *niters);
void setup_polymer_field(double *x,int);
void setup_polymer_simple(double *x, int i);
void setup_polymer_rho(double *x,int);
void setup_polymer_G(double *x);
void setup_const_density(double *x, double *rho,int nloop,int index);
void setup_exp_density(double *x, double *rho,int nloop,int index);
void setup_step_2consts(double *x);
void setup_linear_profile(double *x);
void setup_rho_bar(double *x);
void setup_elec_pot(double *x,int iguess);
void setup_chem_pot(double *x);
int find_length_of_file(char *filename);
void read_in_a_file(int iguess,char *filename);
void shift_the_profile(double *x_new,double fac);
int locate_inode_old(int *ijk);
void communicate_profile(double *x_new,double *x);
void check_zero_densities(double *x);
void chop_profile(double *x, int iguess);
double interpolate(double *pos,int icomp, int nodes_old,
		   double **pos_old, double *Esize_old);
