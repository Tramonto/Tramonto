/*
  Prototypes for dft_vext.c
*/
void setup_external_field_n (int **nelems_w_per_w, int ***elems_w_per_w);
void read_external_field_n (char *filename);
void read_zero_density_TF (char *filename);
void setup_zero ();
void setup_vext_max ();
void setup_semiperm (int **nelems_w_per_w, int ***elems_w_per_w);
void setup_1Dvext (int);
void setup_1Dvext_xmin (int);
void setup_vext_exp (int);
void setup_vext_LJ_atomic (int);
void setup_vext_HS_atomic (int);
void setup_wall_wall (int **nelems_w_per_w, int ***elems_w_per_w,
		      double **uww_tmp);
void setup_integrated_LJ_walls (int, int *, int **);
void comm_wall_els (int, int **nelems_w_per_w, int ***elems_w_per_w,
		    int *nelems_w_per_w_global, int **elems_w_per_w_global);
#ifdef PARALLEL
void correct_zeroTF_array ();
#endif
void comm_vext_max (int *nnodes_vext_max, int **nodes_vext_max);
void comm_loc_to_glob_vec (int *n_loc, int *in_loc_vec, int *out_glob_vec);
