/*
  Prototypes for functions in dft_comm.h
*/
void gather_global_vec (double *loc_vec, int *loc_index, int N_loc,
			double *global_vec);
void gather_global_vec_int (int *loc_vec, int *loc_index, int N_loc,
			    int *global_vec);
void share_global_int_vec (int vec_length, int *int_vec);
void share_global_double_vec (int vec_length, double *double_vec);
