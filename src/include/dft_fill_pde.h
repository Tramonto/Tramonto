/*
  Prototypes from dft_fill_pde.c
*/


void set_fem_weights(double **wt_laplace_ptr, double **wt_source_ptr);
void set_fem_1el_weights(double **wt_lp_1el_ptr, double **wt_s_1el_ptr,
                         int ***elem_permute);
void load_poissons_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
                       double *mat_row,
                       double *resid, double *x, int *bindx_tmp, int fill_flag);
void load_polarize_poissons_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
                       double *mat_row,
                       double *resid, double *x, int *bindx_tmp, int fill_flag);
void load_poisson_bc(double *resid,int inode_box,int loc_inode,int loc_i);
void basis_fn_calc(double **phi, double ***grad_phi, double *evol);
void load_nonlinear_transport_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
				  double *mat_row, double *resid, double *x,
				  int *bindx_tmp, int fill_flag, int iunk);
void load_linear_transport_eqn(int i_box, int inode_box, int loc_i, int *ijk_box,
			       double *mat_row, double *resid, double *x,
			       int *bindx_tmp, int fill_flag, int iunk);
