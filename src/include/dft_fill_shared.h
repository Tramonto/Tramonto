/*
  Prototypes from dft_fill_shared.c
*/

void pre_calc_rho_bar (struct RB_Struct *rho_bar, double *x,
		       int fill_flag, int iter,
		       double ***jac_weights_hs, int ***jac_columns_hs,
		       double *fill_time);

double HW_boundary_weight (int icomp, int ilist, double *hw_weight,
			   int inode_box, int *reflect_flag);
