/*
  Prototypes from dft_continuation.c
*/
int solve_continuation (double *x, double *x2, int polymer_flag);
int nonlinear_solver_conwrap (double *x, void *con_ptr, int step_num);
int linear_solver_conwrap (double *x, int jac_flag, double *tmp);
int komplex_linear_solver_conwrap (double *c, double *d,
				   int jac_flag, double *omega, double *tmp);
void mass_matrix_fill_conwrap (double *x, double *rhs);
void mass_matvec_mult_conwrap (double *x, double *y);
void create_shifted_matrix_conwrap ();
void shifted_matrix_fill_conwrap (double sigma);
void shifted_linear_solver_conwrap (double *x, double *y,
				    int jac_flag, double tol);
void destroy_shifted_matrix_conwrap ();
void matrix_residual_fill_conwrap (double *x, double *rhs, int matflag);
void matvec_mult_conwrap (double *x, double *y);
void assign_parameter_tramonto (int cont_type, double param);
void assign_parameter_conwrap (double param);
void assign_bif_parameter_conwrap (double tp_param);
void calc_scale_vec_conwrap (double *x, double *scale_vec, int numUnks);
double gsum_double_conwrap (double sum);
int gmax_int_conwrap (int max);
void random_vector_conwrap (double *x, int numOwnedUnks);
void perturb_solution_conwrap (double *x, double *x_old,
			       double *scale_vec, int numOwnedUnks);
struct con_struct;
void solution_output_conwrap (int num_soln_flag, double *x, double param,
			      double *x2, double param2,
			      double *x3, double param3,
			      int step_num, int num_its,
			      struct con_struct *con);
double free_energy_diff_conwrap (double *x, double *x2);
