/* 
   Prototypes from dft_fill_rosen_rb.c
*/
void load_nonlocal_hs_rosen_rb (int sten_type, int loc_i, int icomp,
				int izone, int *ijk_box, int fill_flag,
				double *x,
				struct RB_Struct *dphi_drb,
				struct RB_Struct *dphi_drb_bulk,
				struct RB_Struct *dphi_drb_bulk_left,
				struct RB_Struct *dphi_drb_bulk_right,
				double *resid, double *mat_row,
				int *bindx_tmp, struct RB_Struct *rho_bar,
				int loc_inode, int resid_only_flag);

void pre_calc_dphi_drb_rb1 (struct RB_Struct *dphi_drb,
			    double *x,
			    struct RB_Struct *dphi_drb_bulk_ptr,
			    struct RB_Struct *dphi_drb_bulk_left_ptr,
			    struct RB_Struct *dphi_drb_bulk_right_ptr,
			    struct RB_Struct *rho_bar, double *fill_time);

void pre_calc_dphi_drb_rb2 (struct RB_Struct *dphi_drb,
			    double *x,
			    struct RB_Struct *dphi_drb_bulk_ptr,
			    struct RB_Struct *dphi_drb_bulk_left_ptr,
			    struct RB_Struct *dphi_drb_bulk_right_ptr,
			    struct RB_Struct *rho_bar, double *fill_time);

void load_rho_bar_s (int sten_type, double *x, int loc_i, int iunk,
		     int loc_inode, int izone, int *ijk_box,
		     int fill_flag, double *resid,
		     double *mat_row, int *bindx_tmp, int resid_only_flag);

void load_rho_bar_v (double *x, int loc_i, int iunk, int loc_inode,
		     int izone, int *ijk_box, int fill_flag,
		     double *resid, double *mat_row, int *bindx_tmp,
		     int resid_only_flag);
