/*
  Prototypes for dft_out_profiles.c
*/
void collect_x_old (double *x, int flag);
void collect_vext_old ();
void print_profile (char *output_file4);
void print_gofr (char *output_file6);
void print_zeroTF (int **zero_TF, char *output_file);
void print_Nodes_to_zone (int *node_to_zone, char *output_file);
void print_charge_surf (double **charge_w_sum, char *output_file);
void print_charge_vol (double *charge_els, char *output_file);
void print_vext (double **vext, char *output_file);
void print_rho_bar (struct RB_Struct *rho_bar, char *output_file);
void print_time_histogram (int *hist, int *niters);
