/*====================================================================
 * ------------------------
 * | CVS File Information |
 * ------------------------
 *
 * $RCSfile$
 *
 * $Author$
 *
 * $Date$
 *
 * $Revision$
 *
 * $Name$
 *====================================================================*/
#ifndef _CON_CONST_H_
#define _CON_CONST_H_

#ifndef lint
static char *cvs_conconsth_id =
  "$Id$";
#endif

/*****************************************************************************/
/*                       DEFINE STATEMENTS                                   */
/*****************************************************************************/
#define TRUE  1
#define FALSE 0

/* Choices for continuation method */
#define ZERO_ORDER_CONTINUATION        0
#define FIRST_ORDER_CONTINUATION       1
#define ARC_LENGTH_CONTINUATION        2
#define TURNING_POINT_CONTINUATION     3
#define PITCHFORK_CONTINUATION         4
#define HOPF_CONTINUATION              5
#define PHASE_TRANSITION_CONTINUATION  6
#define AUGMENTING_CONDITION           7

/* Choices for matrix fill -- what quantities to calculate */
#define  RHS_ONLY        100
#define  MATRIX_ONLY     101
#define  RHS_MATRIX      102

/* Choices for linear solve about the state of the Jacobian matrix*/
#define  NEW_JACOBIAN                200
#define  OLD_JACOBIAN                201
#define  SAME_BUT_UNSCALED_JACOBIAN  202

/* Internal flags for rhs calculation for continuation linear solves*/
#define  CONT_TANGENT      300
#define  ARC_CONT_SOL2     301
#define  TP_CONT_SOL2      302
#define  TP_CONT_SOL3      303
#define  TP_CONT_SOL4      304

/*****************************************************************************/
/*                       STRUCTURE DEFINITIONS                               */
/*****************************************************************************/

/*
 * set_con_struct: structure of input parameters for controllong
 *                 the continuation routines
 */

struct set_con_struct {
   int max_steps; /* Maximum # of continuation steps */
   int max_newton_its; /* Max # of Newton steps, used for step control only */
   double step_ctrl;   /* step aggressiveness -- 0.0 for constant step*/
   double max_param;   /* parameter value to stop at */
   double first_step;  /* Initial step size */
   double tp_param_guess; /* initial value of turning point parameter */
}; 

struct con_struct {
   int method;    /* CONTINUATION or TP_CONTINUATION */
   int numUnks;      /* Number of unknowns on this Proc (including externals)*/
   int numOwnedUnks; /* Number of unknowns on this Proc (NOT incl. externals)*/
   int printproc;    /* Logical indicating if this Proc does print statements*/
   int step_num;     /* Continuation step num; deactivates arclength when 0  */
   double param;     /* value of continuation parameter                      */
   double param_old; /* old value of continuation parameter                  */
   double step;      /* step size of continuation parameter                  */
   double step_old;  /* old step size of continuation parameter              */
   double arc_step;  /* step size of arc length variable                     */
   double delta_param; /* perturbation of parameter for numerical derivs     */
   double dp_ds;     /* derivative of parameter w.r.t. arc length            */
   double tp_param;  /* value of turning point parameter                     */
   double *x;        /* current solution vector                              */
   double *x_old;    /* previous solution vector                             */
   double *x_tang;   /* tangent to the solution vector w.r.t. parameter      */
   double *scale_vec;/* tangent to the solution vector w.r.t. parameter      */
};
/*****************************************************************************/

extern int con_lib(struct con_struct *con, struct set_con_struct *set_con);
extern int arc_length_bordering_alg(double *x, double *delta_x,
		             struct con_struct *con,
		             double reltol, double abstol);
extern int turning_point_alg(double *x, double *delta_x,
		             struct con_struct *con,
			     double reltol, double abstol);
extern int  nonlinear_solver_conwrap(double *x, void *con,
		                     int step_num);
extern void assign_parameter_conwrap(double param);
extern void assign_tp_parameter_conwrap(double tp_param);
extern int  linear_solver_conwrap(double *x, int jac_flag, double *tmp);
extern void calc_scale_vec_conwrap(double *x, double *scale_vec, int numUnks);
extern void exchange_bdry_conwrap(double *x);
extern double gsum_double_conwrap(double sum);
extern void perturb_solution_conwrap(double *x, double *p, double *s, int n);
extern void solution_output_conwrap(double *x, double param, int type,
                                    int stp, int nits);
extern void calc_rhs_continuation(int rhs_type, double *x, double *x_dot,
		           double *scale_vec, double *x_tmp, double dc_p,
                           double con_par, double *r_vec,
			   int num_total_unknowns, int num_owned_unks);
extern void matrix_residual_fill_conwrap(double *x, double *rhs, int matflag);
extern void matvec_mult_conwrap(double *x, double *y);

#endif
