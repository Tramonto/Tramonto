/*
  Prototypes from local_con_bord.c
*/
int arc_length_bordering_alg(double *x, double *delta_x,
		             struct con_struct *con,
		             double reltol, double abstol);
int turning_point_alg(double *x, double *delta_x, struct con_struct *con,
		      double reltol, double abstol);
int pitchfork_alg(double *x, double *delta_x, struct con_struct *con,
	          double reltol, double abstol);
int hopf_alg(double *x, double *delta_x, struct con_struct *con,
	     double reltol, double abstol);
int phase_transition_alg(double *x, double *delta_x,
		         struct con_struct *con, double reltol, double abstol);
void calc_rhs_continuation(int rhs_type, double *x, double *a,
                           double *x_dot, double *scale_vec, double *x_tmp,
                           double con_par, double perturb, double *r_vec,
                           int num_total_unknowns, int num_owned_unks);
int continuation_hook(double *x, double *delta_x, void *con_void,
                      double reltol, double abstol);

