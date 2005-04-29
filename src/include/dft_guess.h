/*
  Prototypes from dft_guess.c
*/

static void set_initial_guess (int iguess,double **xOwned);
static void setup_polymer_field(double **xOwned,int);
static void setup_polymer_simple(double **xOwned, int i);
static void setup_polymer_rho(double **xOwned,int);
static void setup_polymer_G(double **xOwned);
static void setup_const_density(double **xOwned, double *rho,int nloop,int index);
static void setup_exp_density(double **xOwned, double *rho,int nloop,int index);
static void setup_step_2consts(double **xOwned);
static void setup_linear_profile(double **xOwned);
static void setup_rho_bar(double **xOwned);
static void setup_elec_pot(double **xOwned,int iguess);
static void setup_chem_pot(double **xOwned);
static int find_length_of_file(char *filename);
static void read_in_a_file(int iguess,char *filename);
static void shift_the_profile(double *x_new,double fac);
static int locate_inode_old(int *ijk);
static void communicate_profile(double *x_new,double **xOwned);
static void check_zero_densities(double **xOwned);
static void chop_profile(double **xOwned, int iguess);
static double interpolate(double *pos,int icomp, int nodes_old,
                          double **pos_old, double *Esize_old);
