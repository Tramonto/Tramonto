/*
  Prototypes from dft_fill_msr.c
*/
void put_row_in_msr(int i_box, int loc_i,int *bindx, 
                    int *bindx_tmp, int **bindx_2d,
                    double *val, double *mat_row, 
                    int fill_flag,double *x);
void put_zero_in_msr(int loc_i,int **bindx_2d);
void put_euler_lag_in_msr(int i_box, int loc_i, int **bindx_2d);
void put_coarse_in_msr(int mesh_coasen_flag, int i_box, 
		       int loc_i, int **bindx_2d);
void put_1Dsolution_in_msr(int mesh_coasen_flag, int i_box, 
			   int loc_i, int jnode_box,int **bindx_2d);
void put_transport_in_msr(int i_box, int loc_i, int **bindx_2d);
void put_poisson_in_msr(int i_box, int loc_i, int **bindx_2d);
