/*
  Prototypes for local_util.c
*/
void initialize_util_routines(int n_o, int n_t);
void init_vec(double *u);
double *alloc_vec();
void free_vec (double **ptr);
void vec_copy(double *dx, double *dy);
double dp(double *x, double *y);
double ip(double *x, double *y);
double ltransnorm(double *x, double *scale_vec);

