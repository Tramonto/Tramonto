/*
  Prototypes for dft_newton.c
*/

int solve_problem(double **x_internal_ptr, double **x2_internal_ptr,
                  char *output_file1,int *niters);

int newton_solver(double *x, double *x2, double *fill_time, void *con_ptr,
		  int max_iter, double *t_solve_max_ptr, void * aux_info);

int update_solution(int iter,double *x, double *delta_x);

void full_to_msr(double **a, int *bindx, double *val, int maxrow, int maxcol);

double dotp(double *x, double *y);


