/*
  Prototypes from local_con_bord.c
*/
int arc_length_bordering_alg (double *x, double *delta_x,
			      struct con_struct *con,
			      double reltol, double abstol);
int turning_point_alg (double *x, double *delta_x, struct con_struct *con,
		       double reltol, double abstol);
int pitchfork_alg (double *x, double *delta_x, struct con_struct *con,
		   double reltol, double abstol);
int hopf_alg (double *x, double *delta_x, struct con_struct *con,
	      double reltol, double abstol);
int phase_transition_alg (double *x, double *delta_x,
			  struct con_struct *con, double reltol,
			  double abstol);
void calc_rhs_continuation (int rhs_type, double *x, double *resid_vector,
			    double *ab_vec, double *scale_vec, double *x_tmp,
			    double param, double *r_vec, int numUnks,
			    int numOwnedUnks);
int continuation_hook (double *x, double *delta_x, void *con_void,
		       double reltol, double abstol);
