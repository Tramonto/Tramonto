/*
  Prototypes from dft_fillp.c
*/

void fill_resid_and_matrix_P (double *x, double *resid,
			      int **bindx_2d, double *fill_time, int fill_flag,
			      int iter, int resid_only_flag,int unk_flag);
double load_polymer_G(int,int,int,int,int,int *,int,double *,int *,double *);
double load_polymer_cr(int sten_type, int loc_i, int itype_mer, 
		       int izone, int *ijk_box, int fill_flag,
		       double *resid, double *mat_row, int *bindx_tmp, double *x);
