/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

/* DOCUMENTATION FOR CONTINUATION LIBRARY:    10/26/1999
 *    Andrew Salinger  Org. 9221, MS-1111, 845-3523
 *
 * The continuation library currently consists of 4 files. 
 * (1) [rf_continuation.c]  The MPSalsa version of the interface
 *     to the continuation library. This file (and only this file)
 *     must be extensively modified for every new code that
 *     uses the library. It consists of solve_continuation, the
 *     top level routine called by the application code, and ~10
 *     wrapper routines, which must be filled for your specific
 *     code. The most noteworthy of these include
 *     nonlinear_solver_conwrap, linear_solver_conwrap,
 *     matrix_residual_fill_conwrap, and assign_parameter_conwrap.
 *     Each wrapper has its own comments.
 * (2) [rf_con_lib.c]  This is the stepping algorithm for zeroth-order,
 *     first-order, and arc-length continuation. The turning point
 *     tracking algorithm is only zeroth-order. This routine includes
 *     such things as step-size control and predictor steps. Ideally,
 *     this file does not need to be modified when linking to a new
 *     application code.
 * (3) [rf_con_bord.c] This file contains the bordering algorithms,
 *     currently for arc-length continuation and turning point tracking.
 *     These routines are accessed from within the Newton iteration,
 *     and make frequent calls to the wrapper routines (matrix fill
 *     and solve). Ideally, this file does not need to be modified when
 *     linking to a new application code.
 * (4) [rf_con_const.h] This header file includes definitions of the
 *     continuation structures, define statements for flags used
 *     by the continuation library, and prototypes for the wrapper
 *     routines. Ideally, this file does not need to be modified when
 *     linking to a new application code.
 *
 * How to interface to this library:
 * (0) Have a steady-state code that uses Newton's method
 * (1) Call solve_continuation from your code, in the place
 *     where you normally call the steady-state or transient
 *     drivers.
 * (2) In your current nonlinear solver, add "(void *) con_ptr"
 *     to the argument list. Pass NULL in that place for all
 *     calls to the nonlinear solver besides the one from
 *     nonlinear_solver_conwrap below. In your Newton loop,
 *     after the matrix equation has been solved for the
 *     update vector, delta_x, but before the solution vector, x,
 *     has been updated by delta_x, put in the following call:
 *      if (con_ptr != NULL) continuation_converged =
 *         continuation_hook(x, delta_x, con_ptr, Reltol, Abstol);
 *     (Reltol=1.0e-3, Abstol=1.0e-8 are good defaults.)
 *     When checking convergence of your Newton's method, also
 *     check that (continuation_converged == TRUE).
 * (3) Change the contents of all routines in this file. In
 *     solve_continuation, many fields of the con and set_con structures
 *     must be set. Follow the template and the comments.
 *     Also, the passdown structure can be defined and set with 
 *     any information that needs to be passed from the 
 *     solve_continuation routine to the wrapper routines below
 *     that isn't needed by the continuation library. All
 *     of the wrapper routines (all routines in this file besides
 *     solve_continuation, which all have names ending in _conwrap)
 *     must be filled in with the corresponding call to your
 *     application code. Each has its own comments.
 */

/* Put include statements for your code here. */

/*****************************************************************************/

/* This include file is for the continuation structures, defines, prototypes */

#include <stdio.h>
#include "loca_const.h"
#include "dft_continuation.h"

static void print_con_struct(const struct con_struct* con);

/* Define passdown structure: this structure is global to this file, and
 * provides a way to pass variables from the top solve_continuation routine
 * down to the various wrapper routines below, without passing through the
 * continuation library. For instance, the Jacobian matrix is never seen
 * in the continuation routines, but is needed by several of the wrapper
 * routines in this file. The passdown structure is a way make the
 * matrix global to this file.
 */

struct passdown_struct {
double epswf_previous; /* Save old value of Eps_wf for rescaling external field*/
double epsff_previous; /* Save old value of Eps_ff for rescaling stencils */
double temp_previous; /* Save old value of Temp for rescaling external field*/
double chg_scale_previous;
double **xBox;
double **xOwned;
} passdown;

double get_init_param_value(int cont_type,int);
static void print_final(double param, int step_num);
static void translate_2dBox_1dOwned(double **xBox, double *x);
static void translate_1dOwned_2dBox(double *x, double **xBox);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int solve_continuation( double **xx, double **xx2)

/* Interface routine to the continuation library.
 */
{
  struct con_struct con;
  int nstep=0;
  double *x, *x2;


  /******************************* First Executable Statment *****************/

  /* xOwned and xBox are always temp storage space */
  passdown.xOwned = (double **) array_alloc(2, Nunk_per_node, Nnodes_per_proc, sizeof(double));
  passdown.xBox   = (double **) array_alloc(2, Nunk_per_node, Nnodes_box,      sizeof(double));
  x               = (double * ) array_alloc(1, Nunk_per_node*Nnodes_per_proc,  sizeof(double));

  /* Translate from 2d to 1d data structure */
  translate_2dBox_1dOwned(xx, x);

  /* 
   * Set passdown structure -- variables needed in the argument
   * lists to wrapped routines but not needed in the continuation
   * library.
   */

  /*
   * Fill all the necessary information into the 'con' structure
   * for running the continuation routines. (In this example, this
   * info was read into an identical structure, so the copying isn't
   * necessary. However, if a application has input the continuation
   * information in a different format, these lines below will be
   * needed to load the 'con' structure.)
   */

  /* First load the general info structure */

  switch (Loca.method) {
    case 0: con.general_info.method = ZERO_ORDER_CONTINUATION; break;
    case 1: con.general_info.method = FIRST_ORDER_CONTINUATION; break;
    case 2: con.general_info.method = ARC_LENGTH_CONTINUATION; break;
    case 3: con.general_info.method = TURNING_POINT_CONTINUATION; break;
    case 4: con.general_info.method = PHASE_TRANSITION_CONTINUATION; break;
    default:
       if (Proc==0 && Iwrite_screen !=SCREEN_NONE) 
          printf("\nERROR in solve_continuation: Unknown " "continuation method: %d\n",Loca.method);
       exit(-1);
  }
 
  /* Allocate second solution vector for phase transitions */
  if (con.general_info.method==PHASE_TRANSITION_CONTINUATION) {
    x2              = (double * ) array_alloc(1, Nunk_per_node*Nnodes_per_proc,  sizeof(double));
    translate_2dBox_1dOwned(xx2, x2);
  }

  if (Iwrite_screen==SCREEN_VERBOSE) printf("\n###\nPROC = %d\n\n",Proc);

  con.general_info.param        = get_init_param_value(Loca.cont_type1,0);
  con.general_info.x            = x;
  con.general_info.numUnks      = Nunk_per_node*Nnodes_per_proc;
  con.general_info.numOwnedUnks = Nunk_per_node*Nnodes_per_proc;
  if (Proc==0) {
     switch (Iwrite_screen) {
	     case SCREEN_NONE: con.general_info.printproc = 0; break;
	     case SCREEN_ERRORS_ONLY:   con.general_info.printproc = 1; break;
	     /*case SCREEN_BASIC: con.general_info.printproc = 5; break;*/
	     case SCREEN_BASIC: con.general_info.printproc = 1; break;
	     case VERBOSE:   con.general_info.printproc = 8; break;
     }
  }
  else         con.general_info.printproc = FALSE;
  con.general_info.nv_save   = FALSE;   /* Salsa saves null vector anyway */
  con.general_info.perturb  = 1.0e-7;


  /* Then load the stepping info structure */

  con.stepping_info.first_step     = Loca.step_size;
  con.stepping_info.base_step      = 0;
  con.stepping_info.max_steps      = Loca.num_steps;
  con.stepping_info.max_param      = 1.0e10;
  con.stepping_info.max_delta_p    = 1.0e10;
  con.stepping_info.min_delta_p    = 1.0e-10;
  con.stepping_info.step_ctrl      = Loca.aggr;
  con.stepping_info.max_newton_its = Max_NL_iter;

  /* Then load one of the method dependent structures */

  switch (con.general_info.method) {
    case ZERO_ORDER_CONTINUATION:
    case FIRST_ORDER_CONTINUATION:
      break;
    case ARC_LENGTH_CONTINUATION:
      con.arclength_info.dp_ds2_goal = 0.1;
      con.arclength_info.dp_ds_max   = 0.6;
      con.arclength_info.tang_exp    = 2.0;
      con.arclength_info.tang_step_limit = 0.5;
      break;
    case TURNING_POINT_CONTINUATION:
      con.turning_point_info.bif_param = get_init_param_value(Loca.cont_type2,1);
      con.general_info.nv_restart       = FALSE;
      break;
    case PITCHFORK_CONTINUATION:
      con.pitchfork_info.bif_param = get_init_param_value(Loca.cont_type2,1);
      con.pitchfork_info.psi        = NULL;
      break;
    case HOPF_CONTINUATION:
      con.hopf_info.bif_param = get_init_param_value(Loca.cont_type2,1);
      con.hopf_info.omega       = 0.0;
      con.hopf_info.y_vec       = NULL;
      con.hopf_info.z_vec       = NULL;
      con.hopf_info.mass_flag   = 0;
      break;
    case PHASE_TRANSITION_CONTINUATION:
      con.phase_transition_info.bif_param
          = get_init_param_value(Loca.cont_type2,1);
      con.phase_transition_info.x2  = x2;
      break;
    default:
      if (Iwrite_screen != SCREEN_NONE) printf("ERROR: Unknown LOCA input method: %d\n",con.general_info.method);
      exit(-1);
  }

  /* Finally, load the eigensolver structures */

  con.eigen_info.Num_Eigenvalues   = 0;
  con.eigen_info.Shift_Point[0]    = 5.0;
  con.eigen_info.Shift_Point[1]    = 50.0;
  con.eigen_info.Shift_Point[2]    = 1.0;
  con.eigen_info.Arnoldi           = 30;
  con.eigen_info.Residual_Tol[0]   = 1.0e-4;
  con.eigen_info.Residual_Tol[1]   = 1.0e-4;
  con.eigen_info.Max_Iter          = 1;
  con.eigen_info.Every_n_Steps     = 1;

  /* print out continuation structure */

  if (Iwrite_screen == SCREEN_VERBOSE) print_con_struct(&con);

  /* Now call continuation library and return */

  nstep = con_lib(&con, NULL);


  /*****  Clean Up ********/
  /* Load final solutions back into original memory location */
  translate_1dOwned_2dBox(x, xx);
  if (con.general_info.method==PHASE_TRANSITION_CONTINUATION)
      translate_1dOwned_2dBox(x2, xx2);

  /* Free temporary memory */
  safe_free((void **) &passdown.xOwned);
  safe_free((void **) &passdown.xBox);
  safe_free((void **) &x);
  if (con.general_info.method==PHASE_TRANSITION_CONTINUATION)
    safe_free((void **) &x2);

  if (con.general_info.printproc && Iwrite_screen == SCREEN_VERBOSE) print_final(con.general_info.param, nstep);

  return nstep;
} /**************** END of solve_continuation () *****************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int nonlinear_solver_conwrap(double *x, void *con_ptr, int step_num,
		             double lambda, double delta_s, void * aux_info)
/* Put the call to your nonlinear solver here.
 * Input:
 *    x         solution vector
 *    con_ptr   pointer to continuation structure, cast to (void *)
 *              must be passed to nonlinear solver and then passed
 *              to bordering algorithms.
 *    step_num  Continuation step number
 *
 * Output:
 *
 * Return Value:
 *    num_its  Number of Newton iterations needed for
 *             convergence, used to pick next step size.
 *             Negative value means nonlinear solver didn't converge.
 */
{
  int num_its;

  translate_1dOwned_2dBox(x, passdown.xBox);

  num_its = newton_solver(passdown.xBox, con_ptr);

  translate_2dBox_1dOwned(passdown.xBox, x);

  return (num_its);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int linear_solver_conwrap(double *x, int jac_flag, double *tmp)
/* Put the call to your linear solver here. There are three options
 * about the reuse of the preconditioner. It is always safe to 
 * just solve the matrix from scratch.
 * Input:
 *    x          Right hand side
 *    jac_flag   Flag indicating the status of the Jacobian so that
 *               preconditioners can be used: 
 *               NEW_JACOBIAN:   recalculate preconditioner
 *               OLD_JACOBIAN:   reuse preconditioner
 *               SAME_BUT_UNSCALED_JACOBIAN: Must rescale the matrix and
 *                               can reuse preconditioner. This happens
 *                               when the matrix has been recalculated
 *                               at the same conditions as before.
 *    tmp        Work space array same length as x, only used for
 *               the SAME_BUT_UNSCALED_JACOBIAN option.
 *
 * Output:
 *    x          Solution vector
 *
 * Return Value:
 *    Negative value means linear solver didn't converge.
 */
{
  int i,j;

  for (i=0; i<Nunk_per_node; i++)
    for (j=0; j<Nnodes_per_proc; j++)
      passdown.xOwned[i][j] = x[i*Nnodes_per_proc + j];
  (void) dft_linprobmgr_setrhs(LinProbMgr_manager, passdown.xOwned);

  if (jac_flag == OLD_JACOBIAN || jac_flag == CHECK_JACOBIAN) {
    /* reuse the preconditioner for this same matrix */
    /*Aztec.options[AZ_pre_calc] = AZ_reuse;*/
  }
  else if (jac_flag == SAME_BUT_UNSCALED_JACOBIAN) {
  }
  else if (jac_flag != NEW_JACOBIAN) {
    if (Iwrite_screen != SCREEN_NONE) printf("ERROR: linear solve conwrap: unknown flag value %d\n",jac_flag);
    exit(-1);
  }

 (void) dft_linprobmgr_setupsolver(LinProbMgr_manager);
 (void) dft_linprobmgr_solve(LinProbMgr_manager);
 (void) dft_linprobmgr_getlhs(LinProbMgr_manager, passdown.xBox);

 translate_2dBox_1dOwned(passdown.xBox, x);

  return 0;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int komplex_linear_solver_conwrap(double *c, double *d,
    int jac_flag, double *omega, double *tmp){return 0;}
void mass_matrix_fill_conwrap(double *x, double *rhs){}
void mass_matvec_mult_conwrap(double *x, double *y){}
void create_shifted_matrix_conwrap(){}
void shifted_matrix_fill_conwrap(double sigma){}
void shifted_linear_solver_conwrap(double *x, double *y,
                                   int jac_flag, double tol){}
void destroy_shifted_matrix_conwrap(){}

/*****************************************************************************/
void matrix_residual_fill_conwrap(double *x, double *rhs, int matflag)
/* Put the call to your matrix/residual fill routine here.
 * Input:
 *    x         Solution vector
 *    matflag   Flag indicating residual (RHS_ONLY), matrix (MATRIX_ONLY),
 *              or both (RHS_MATRIX) are requested.
 *
 * Output:
 *    rhs       Right hand side
 *
 * Return Value:
 */
{
  int i, j, resid_only_flag;
  double l2_resid,resid_sum;

  if (matflag == RHS_ONLY) {
      resid_only_flag = TRUE;
      (void) dft_linprobmgr_initializeproblemvalues(LinProbMgr_manager);
  }
  else { 
      resid_only_flag = FALSE;
      (void) dft_linprobmgr_initializeproblemvalues(LinProbMgr_manager);
  }


  /*fill_time not currently plugged in, iter hardwire above 2 */
  translate_1dOwned_2dBox(x, passdown.xBox);

  resid_sum=fill_resid_and_matrix_control(passdown.xBox, 0, resid_only_flag);
  (void) dft_linprobmgr_finalizeproblemvalues(LinProbMgr_manager);

  (void) dft_linprobmgr_getrhs(LinProbMgr_manager, passdown.xOwned);
   for (i=0; i<Nunk_per_node; i++)
     for (j=0; j<Nnodes_per_proc; j++)
/* NOTE NEGATIVE SIGN TO FIT CONVENTION !!!! */
       rhs[i*Nnodes_per_proc + j] = -passdown.xOwned[i][j];


  /* calculate and print resid norm */
  l2_resid=0.0;
  for (i=0; i<Nunk_per_node*Nnodes_per_proc; i++) l2_resid += rhs[i]*rhs[i];

  l2_resid= sqrt(gsum_double_conwrap(l2_resid));
  if (Proc==0 && Iwrite_screen ==SCREEN_VERBOSE) printf("\t\tNorm of resid vector = %20.15g\n", l2_resid);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void matvec_mult_conwrap(double *x, double *y)
/* Put the call to your matrix-vector multiply here.
 * Input:
 *    x         Vector of length number of unknowns
 *
 * Output:
 *    y         Matrix times x.
 *
 * Return Value:
 */
{
  int i,j;

  translate_1dOwned_2dBox(x, passdown.xBox);
  (void) dft_linprobmgr_applymatrix(LinProbMgr_manager,passdown.xBox, passdown.xOwned);

   for (i=0; i<Nunk_per_node; i++)
     for (j=0; j<Nnodes_per_proc; j++)
       y[i*Nnodes_per_proc + j] = passdown.xOwned[i][j];

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_parameter_conwrap(double param)
/* Put the call to a routine to assign the continuation parameter here.
 * Input:
 *    param     New value of continuation parameter.
 *
 * Output:
 *
 * Return Value:
 */
{
  if (Proc==0 && Iwrite_screen !=SCREEN_NONE && Iwrite_screen !=SCREEN_ERRORS_ONLY) {
        printf("\tContinuation parameter #%d set to %g\n", Loca.cont_type1, param);
  }
  assign_parameter_tramonto(Loca.cont_type1, param,0);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void assign_bif_parameter_conwrap(double tp_param)
/* Put the call to a routine to assign the continuation parameter here.
 * Input:
 *    tp_param     New value of continuation parameter.
 *
 * Output:
 *
 * Return Value:
 */
{
  if (Proc==0 && Iwrite_screen ==SCREEN_VERBOSE){ 
        printf("\tSecond (floating) parameter #%d set to %20.15g\n", Loca.cont_type2, tp_param);
  }
  assign_parameter_tramonto(Loca.cont_type2, tp_param,1);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void calc_scale_vec_conwrap(double *x, double *scale_vec, int numUnks)
/* Put the call to a routine to calculate a scaling vector here.
 * Input:
 *    x          New value of continuation parameter.
 *    numUnks    Number of unknowns on this proc, the length of x
 *               and scale_vec.
 *
 * Output:
 *    scale_vec  Vector of length number of unknowns used to scale
 *               variables so that one type of unknown (e.g. pressure)
 *               doesn't dominate over others. Used to balance the
 *               variables and the arc-length variable in arc-length
 *               continuation, and for scaling the null vector in
 *               turning point tracking. Using reciprocal of the average
 *               value of that variable type is a good choice. Vector
 *               of all ones should suffice for most problems.
 *
 *               For tramonto, I've seen better behavior with vector 
 *               of all ones instead of using Rho_b.
 *
 * Return Value:
 */
{
  int      i;

  /* This was the default -- changed for spinodal tracking */
  for (i=0; i<numUnks; i++) scale_vec[i] = 1.0;
  /*for (i=0; i<numUnks; i++) scale_vec[i] = 1.0/ (sqrt(fabs(x[i])) + 1.0e-5);*/

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double gsum_double_conwrap(double sum)
/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    sum     Value of double on this processor to be summed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global sum is returned on all processors.
 */
{
  return gsum_double(sum);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double gmax_double_conwrap(double sum)
{ return gmax_double(sum); }

/* This worked calling the C code from C++. Probably not needed. */
void fill_resid_and_matrix_control_conwrap(double** xBox, int ii, int jj)
{  double r;
   r= fill_resid_and_matrix_control(xBox, ii, jj); 
}

void safe_free_conwrap(void** p)
{  safe_free(p); }

double** array_alloc_2d_conwrap(unsigned int ii, unsigned int jj, unsigned int kk)
{  return (double **) array_alloc_2d(ii, jj, kk); }
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int gmax_int_conwrap(int max)
/* Put the call to a routine to calculate a global sum.
 * Just return sum for single processor jobs.
 * Input:
 *    max     Value of integer on this processor to be maxed on all procs.
 *
 * Output:
 *
 * Return Value:
 *    The global max is returned on all processors.
 *
 * Only used by Eigensolver
 */
{
  return gmax_int(max);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void random_vector_conwrap(double *x, int numOwnedUnks)
/* Put a routine to calculate a random vector.
 * Input:
 *    numOwnedUnks  Length of owned nodes part of x.
 *
 * Output:
 *    x             Random vector.
 *
 * Used by eigensolver only
 */
{
  int i;
  if (Iwrite_screen != SCREEN_NONE) printf("\tWARNING: random_vector_conwrap just filling vec with constant\n");
  for (i=0; i<numOwnedUnks; i++) x[i] = 0.5;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void perturb_solution_conwrap(double *x, double *x_old,
		              double *scale_vec, int numOwnedUnks)
/* Put a routine to perturb the solution vector a little bit here.
 * This is to move a converged solution at a singularity off
 * of the singularity before doing continuation. This ain't pretty
 * but has helped convergence on some turning point tracking problems.
 * Input:
 *    x_old         Current solution vector.
 *    scale_vec     Work space for a vector to scale x.
 *    numOwnedUnks  Length of owned nodes part of x, x_old, scale_vec
 *
 * Output:
 *    x             Solution vector perturbed a bit.
 *    
 * Return Value:
 */
{
  int i;

  random_vector_conwrap(x, numOwnedUnks);

  for (i=0; i<numOwnedUnks; i++) x[i] = x_old[i] * (1.0 + x[i]*1.0e-5);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void solution_output_conwrap(int num_soln_flag, double *x, double param,
                             double *x2, double param2,
                             double *x3, double param3,
                             int step_num, int num_its, struct con_struct *con)

/* Put the call to your solution output (both file and screen) routines here.
 * Input:
 *    num_soln_flag  Flag for how many solution vectors are being passed for
 *                   output. For parameter continuation and turning point
 *                   tracking there is just 1, for pitchfork and phase
 *                   transitions there are 2, and for Hopfs there are 3.
 *                   The eigensolver sends -1 for a single eigenvector and
 *                   -2 fir a complex pair of eigenvectors.
 *    x            First solution vector for output (x or first eigen.
 *    param        Continuation parameter value (Real part ev when flag < 0)
 *    x2           Second solution vector for output (y_vec or x2)
 *    param2       Bifurcation parameter value (imaginary part ev if flag=-2)
 *    x3           Third solution vector for output (z_vec for Hopfs)
 *    param3       Third Parameter value (frequency Hopfs)
 *    step_num+1   Time index to output to (step_num is 0 based).
 *    num_its      Number of Newton iterations used for for convergence
 *    con          pointer to continuation structure, for passing to
 *                 the eigensolver.
 *
 * Output:
 *
 * Return Value:
 */
{
  double time_save=0;

  translate_1dOwned_2dBox(x, passdown.xBox);

  post_process(passdown.xBox, &num_its, &time_save,step_num, FALSE,FROM_LOCA);

  /* Print second solution for Phase Trans, but not Spinodal tracking */

  if (num_soln_flag == 2 
     && con->general_info.method == PHASE_TRANSITION_CONTINUATION) {

    translate_1dOwned_2dBox(x2, passdown.xBox);

    post_process(passdown.xBox, &num_its, &time_save, step_num, TRUE,FROM_LOCA);
  }
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/* Get around some C++ vs C linkage issue */
double calc_free_energy_conwrap(double **xB)
{ return calc_free_energy(NULL, xB); }


double free_energy_diff_conwrap(double *x, double *x2)
/* Call to return the free energy difference betwen two solutions
 * Input:
 *    x    One solution vector
 *    x2   Second solution vector
 *
 * Output:
 *
 * Return Value:
 *    The difference in the free energy beween the two solutions
 */
{
  double energy1, energy2;

  setup_integrals();

  translate_1dOwned_2dBox(x, passdown.xBox);
  energy1 = calc_free_energy(NULL,passdown.xBox);

  translate_1dOwned_2dBox(x2, passdown.xBox);
  energy2 = calc_free_energy(NULL,passdown.xBox);

  if (Iwrite_screen == SCREEN_VERBOSE) printf("energy1 %12.8g energy2  %12.8g\n",energy1, energy2);
  return energy1-energy2;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_con_struct(const struct con_struct* con)
/* Routine for printing out the con structure to the screen */

{
  printf("--------------------------------------------------\n"); 
  printf("\tcon->general_info.param=             %10.4g\n",
    con->general_info.param);
  printf("\tcon->general_info.numUnks=           %10d\n",
    con->general_info.numUnks);
  printf("\tcon->general_info.numOwnedUnks=      %10d\n",
    con->general_info.numOwnedUnks);
  printf("\tcon->general_info.printproc=         %10d\n",
    con->general_info.printproc);

  printf("\tcon->stepping_info.first_step=       %10.4g\n",
    con->stepping_info.first_step);
  printf("\tcon->stepping_info.max_steps=        %10d\n",
    con->stepping_info.max_steps);
  printf("\tcon->stepping_info.max_param=        %10.4g\n",
    con->stepping_info.max_param);
  printf("\tcon->stepping_info.max_delta_p=      %10.4g\n",
    con->stepping_info.max_delta_p);
  printf("\tcon->stepping_info.step_ctrl=        %10.4g\n",
    con->stepping_info.step_ctrl);
  printf("\tcon->stepping_info.max_newton_its=   %10d\n",
    con->stepping_info.max_newton_its);

  if(con->general_info.method==ARC_LENGTH_CONTINUATION) {
    printf("\tcon->arclength_info.dp_ds2_goal=     %10.4g\n",
      con->arclength_info.dp_ds2_goal);
    printf("\tcon->arclength_info.dp_ds_max=       %10.4g\n",
      con->arclength_info.dp_ds_max);
    printf("\tcon->arclength_info.tang_exp=        %10.4g\n",
      con->arclength_info.tang_exp);
    printf("\tcon->arclength_info.tang_step_limit= %10.4g\n",
      con->arclength_info.tang_step_limit);
  }

  if(con->general_info.method==TURNING_POINT_CONTINUATION)
    printf("\tcon->turning_point_info.bif_param=   %10.4g\n",
      con->turning_point_info.bif_param);

  if(con->general_info.method==PITCHFORK_CONTINUATION)
    printf("\tcon->pitchfork_info.bif_param=       %10.4g\n",
      con->pitchfork_info.bif_param);

  if(con->general_info.method==HOPF_CONTINUATION)
    printf("\tcon->hopf_info.bif_param=            %10.4g\n",
      con->hopf_info.bif_param);

  if(con->eigen_info.Num_Eigenvalues > 0) {
    printf("\tcon->eigen_info.Shift_Point[0]=      %10.4g\n",
      con->eigen_info.Shift_Point[0]);
    printf("\tcon->eigen_info.Shift_Point[1]=      %10.4g\n",
      con->eigen_info.Shift_Point[1]);
    printf("\tcon->eigen_info.Arnoldi=             %10d\n",
      con->eigen_info.Arnoldi);
    printf("\tcon->eigen_info.Residual_Tol[0]=     %10.4g\n",
      con->eigen_info.Residual_Tol[0]);
    printf("\tcon->eigen_info.Every_n_Steps=       %10d\n",
      con->eigen_info.Every_n_Steps);
  }
  printf("--------------------------------------------------\n"); 
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_final(double param, int step_num)

/*
 * Print out the final results and counters
 */

{
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  printf("CONTINUATION ROUTINE HAS FINISHED: \n");
  printf("\tEnding Parameter value     = %g\n", param);
  printf("\tNnumber of steps           = %d\n", step_num);
  printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
} 
/************* END of print_final () ***************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void translate_1dOwned_2dBox(double *x, double **xBox) 
{
  int i,j;
  for (i=0; i<Nunk_per_node; i++)
    for (j=0; j<Nnodes_per_proc;j++)
      passdown.xOwned[i][j] = x[i*Nnodes_per_proc + j];
  (void) dft_linprobmgr_importr2c(LinProbMgr_manager, passdown.xOwned, xBox);
}

static void translate_2dBox_1dOwned(double **xBox, double *x) 
{
  int i,j;
  box2owned(xBox, passdown.xOwned);
  for (i=0; i<Nunk_per_node; i++)
    for (j=0; j<Nnodes_per_proc; j++)
      x[i*Nnodes_per_proc + j] = passdown.xOwned[i][j];
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int continuation_hook_conwrap(double **xx, double **delta_xx, void *con_ptr,
		              double reltol, double abstol)
{
  int converged,i;
  double *x, *delta_x;
  x       = (double *) array_alloc(1, Nunk_per_node*Nnodes_per_proc,  sizeof(double));
  delta_x = (double *) array_alloc(1, Nunk_per_node*Nnodes_per_proc,  sizeof(double));

  translate_2dBox_1dOwned(xx, x);
  translate_2dBox_1dOwned(delta_xx, delta_x);

    for (i=0; i<Nunk_per_node*Nnodes_per_proc; i++) delta_x[i] *= -1.0;
    converged = continuation_hook(x, delta_x, con_ptr, reltol, abstol);
    for (i=0; i<Nunk_per_node*Nnodes_per_proc; i++) delta_x[i] *= -1.0;

  translate_1dOwned_2dBox(x, xx);
  translate_1dOwned_2dBox(delta_x, delta_xx);

  safe_free((void **) &x);
  safe_free((void **) &delta_x);

  return converged;
}
