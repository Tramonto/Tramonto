/*
  Prototypes from dft_fill_rosen.c
*/

double load_nonlocal_hs_rosen (int sten_type, int loc_i, int icomp,
			       int izone, int *ijk_box, int fill_flag,
			       struct RB_Struct *rho_bar,
			       struct RB_Struct *dphi_drb,
			       struct RB_Struct *dphi_drb_bulk,
			       struct RB_Struct *dphi_drb_bulk_left,
			       struct RB_Struct *dphi_drb_bulk_right,
			       double *resid, double *mat_row, int *bindx_tmp,
			       double ***jac_weights_hs,
			       int ***jac_columns_hs, int resid_only_flag);

void pre_calc_dphi_drb (struct RB_Struct *dphi_drb,
			struct RB_Struct *rho_bar,
			struct RB_Struct *dphi_drb_bulk_ptr,
			struct RB_Struct *dphi_drb_bulk_left_ptr,
			struct RB_Struct *dphi_drb_bulk_right_ptr);

void pre_calc_dphi_drb2 (struct RB_Struct *dphi_drb,
			 struct RB_Struct *rho_bar,
			 struct RB_Struct *dphi_drb_bulk_ptr,
			 struct RB_Struct *dphi_drb_bulk_left_ptr,
			 struct RB_Struct *dphi_drb_bulk_right_ptr);

struct RB_Struct d2phi_drb2_delta (struct RB_Struct rho_bar, double weight,
				   int *offset, double inv_rad,
				   double inv_4pir, double inv_4pirsq,
				   double *sign, int fill_flag);

struct RB_Struct d2phi_drb2_delta2 (struct RB_Struct rho_bar, double weight,
				    int *offset, double inv_rad,
				    double inv_4pir, double inv_4pirsq,
				    double *sign, int fill_flag);

struct RB_Struct d2phi_drb2_theta (struct RB_Struct rho_bar, double weight,
				   int *offset, int fill_flag);

struct RB_Struct d2phi_drb2_theta2 (struct RB_Struct rho_bar, double weight,
				    int *offset, int fill_flag);

void load_Jac_Save_Jacobian (double *mat_row, struct RB_Struct tmp,
			     int jnode_box, int jzone, int *bindx_tmp,
			     int fill_flag, double ***jac_weights_hs,
			     int ***jac_columns_hs);

void load_Jacobian (double *mat_row, struct RB_Struct tmp, int jnode_box,
		    int jzone, int *bindx_tmp, int fast_fill_TF_j,
		    int fill_flag);

void fast_fill_Jacobian (int sten_type,
			 int jcomp, int jzone, int jnode_box,
			 double tmp_scalar, double *tmp_vector,
			 double *mat_row, int *bindx_tmp, int fill_flag);

void medium_fill_Jacobian (int sten_type,
			   int jcomp, int jzone, int jlist, int jnode_box,
			   double tmp_scalar, double *tmp_vector,
			   double *mat_row, int *bindx_tmp, int fill_flag);

void slow_fill_Jacobian (int sten_type,
			 int jcomp, int jzone, int jlist, int jnode_box,
			 double tmp_scalar, double *tmp_vector,
			 double *mat_row, int *bindx_tmp, int fill_flag);
