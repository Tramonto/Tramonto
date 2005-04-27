/*
  Prototypes for dft_newton.c
*/

int solve_problem(double **x_internal_ptr, double **x2_internal_ptr);

int newton_solver(double **x, void *con_ptr);

int update_solution(double **x, double **delta_x, int iter);

/* COnversion routines */
void internal2local(double *xInternal, double** xLocal);
void local2internal(double** xLocal, double *xInternal);
void box2local(double** xBox, double** xLocal);

