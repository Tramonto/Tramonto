/*
  Prototypes for dft_newton.c
*/

int solve_problem(double **x, double **x2);

int newton_solver(double **x, void *con_ptr);

int update_solution(double **x, double **delta_x, int iter);

/* Conversion routine */
void box2owned(double** xBox, double** xOwned);

