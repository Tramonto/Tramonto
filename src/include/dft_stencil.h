/*
  Prototypes for dft_stencil.c
*/

void calc_stencils(void);
void shorten_stencil(struct Stencil_Struct *sten);
int ijk_to_isten_index(int *ijk,int *el_in_radius);
void sort2_int_double_array(int n, int ra[], double rb[],
				   int *rc[], int ndim);
void set_gauss_quad(int ngp, double *gp, double *gw);
double int_cr(double r_low,double r_upp,double slope_dr,int icomp,int jcomp,
	      int irmin, double zsq, double *rx_low);
double gauss(double r, int i, int j);
void gser(double *gamser, double a, double x, double *gln);
void gcf(double *gammcf, double a, double x, double *gln);
double gammln(double xx);

