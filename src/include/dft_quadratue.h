/*
  Prototypes from dft_quadrature.c
*/
int get_integration_pts (int isten, int izone,
			 double ***point_ptr, double **wt_ptr);
void delta_tetrahedron (double **point, double *wt);
void delta_octahedron (double **point, double *wt);
void delta_cube (double **point, double *wt);
void delta_icosahedron (double **point, double *wt);
void delta_dodecahedron (double **point, double *wt);
void delta_one_D_twelve (double **point, double *wt);
void delta_midpoint (double **point, double *wt, int izone);
void theta_midpoint (double **point, double *wt, int izone, int num_dim);
void delta_seventy_two (double **point, double *wt);
void get_radial_quadrature (double gauss_pt[], double gauss_wt[], int num_gp);
