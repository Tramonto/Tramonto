/*====================================================================
 ------------------------
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
#ifndef lint
static char *cvs_conlib_id =
  "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
#include "rf_salsa.h"
*/

#include "rf_allo.h"
#include "rf_con_const.h"

/******************************************************************************
*
*           GLOBALS DEFINED IN THIS FILE
*
*       Function                type                       Called By
*       --------------        ---------                ----------------------
*
*    solve_continuation         int                  solve_problem:rf_solve.c
*
******************************************************************************/

/*****************************************************************************/
/********************* Static Functions Used *********************************/
/*****************************************************************************/

double scaled_dot_prod(double *x, double *y, double *sc, int n);
static void print_cont_step1(int order, struct con_struct *con);
static void print_cont_step2(int order, struct con_struct *con);
static void print_cont_step_fail(int order, struct con_struct *con);
static void print_line (char *charstr, int ntimes);
static double simple_step_control(int num_newt_conv, int max_Newton_steps,
  	              		  double step_ctrl);
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int con_lib(struct con_struct *con, struct set_con_struct *set_con)

/*****************************************************************************
*
*  Input Variables:
*
******************************************************************************/

{

  /* Local Variables */

  int     n;                    /* Loop index                          */
  int     order;                /* Continuation order flag: 
				  	  0 - zero-order continuation
					  1 - first-order continuation
					  2 - arc-length continuation
				        This flag is always 0 on the first
					solve, and 0 for turning point or any
					other special continuation */
  int     one = 1, i;
  int     num_newt_conv = 0;    /* Number of newton iterations to reach
                                        convergence for last nonlinear solve
                                        -- used to pick next step size
				        ALSO error flag, when < 0            */
  char    *yo = "con_lib";
  int     sn_old, sn_new;                /* Sign of con->dp_ds, to check for a
                                        turning point                        */
  /******************************* First Executable Statment *****************/

  /*
   * Initialize arrays to store predicted and old solutions. Save the current
   * solution, con->x, into con->x_old
   */

  if (con->method != PHASE_TRANSITION_CONTINUATION &&
      con->method != PITCHFORK_CONTINUATION )
    con->x_tang    = array_alloc (1, con->numUnks, sizeof(double));
  else if (con->x_tang == NULL) {
    printf("ERROR: Second solution vector required for Phase Transition alg"
           "and pitchfork algs\n");
    exit(-1);
  }
  con->x_old     = array_alloc (1, con->numUnks, sizeof(double));
  con->scale_vec = array_alloc (1, con->numUnks, sizeof(double));
  dcopy_(&con->numUnks, con->x, &one, con->x_old, &one);
  con->step      = 0.0;

  /* Adjust the BCs/Properties/whatever that the con param really represents */

  assign_parameter_conwrap(con->param);
  if (con->method == TURNING_POINT_CONTINUATION ||
      con->method == PITCHFORK_CONTINUATION ||
      con->method == PHASE_TRANSITION_CONTINUATION) {
     con->tp_param = set_con->tp_param_guess;
     assign_tp_parameter_conwrap(con->tp_param);
  }

  /* In tp_continuation, perturb initial guess off of potential singularity*/

  if (con->method == TURNING_POINT_CONTINUATION ||
      con->method == PITCHFORK_CONTINUATION) {
    if (con->printproc) printf("\tcon_lib: Adding random"
                               " perturbation for tp_continuation\n");

    perturb_solution_conwrap(con->x, con->x_old,
		             con->scale_vec, con->numOwnedUnks);
  }

  /*
   * Print out general time integration information
   */

  if (con->printproc) {
    printf("\n"); print_line("~", 80); print_line("~", 80);
    printf("%s: Start Continuation\n", yo);
    printf("\tInitial step size = %e \n",set_con->first_step);
    printf("\tMax number of continuation steps = %d\n",
           set_con->max_steps);
    printf("\tMax parameter value = %g\n", set_con->max_param);
    print_line("~", 80); print_line("~", 80); printf("\n");
  }

  /***************************** CONTINUATION LOOP **************************/

  /*
   * Initialize the time step counter to 0. Set order flag to zero-order
   * continuation through first solution.
   */

  con->step_num = order = 0;

  /*
   * Loop through a number of continuation steps - note the loop index may not
   * represent the actual step counter, due to failed steps.  The convention
   * here is that solution 0 is not a step, so there will be
   * set_con->max_steps+1 solutions, numbered 0 through set_con->max_steps.
   */

  for (n = 0; n <= set_con->max_steps; n++) {

    /*
     * Print out an initial statement about the step.
     */

    if (con->printproc) print_cont_step1(order, con);

    /*
    *  Exchange information about con->x
    */

    exchange_bdry_conwrap(con->x);

    /*
     * Solve the system of equations at the current step.
     * Note - con->x is considered to be updated, on return from this
     * solution.
     */

    num_newt_conv = nonlinear_solver_conwrap(con->x,(void *)con,con->step_num);

    /*
     * Check for convergence
     */

    if (num_newt_conv < 0) {

      /*
       * Convergence Failure!
       *
       * If initial guess did not converge, abort entire continuation run
       */

      if (con->step_num == 0) {
        n = set_con->max_steps;       /* Force IO and exit  */
        con->step_num = -1;           /* Set failure flag */

        if (con->printproc) {
          printf("\n\t %s: INITIAL GUESS DID NOT CONVERGE:", yo);
          printf("\n\t\t\t     ABORTING CONTINUATION RUN\n");
        }
      }

      /*
       * If this convergence failure wasn't the first or last step, cut step
       * size in half and calculate a new initial guess.  New guess is
       * the old solution plus the tangent to the previous prediction
       * times the halved step size.
       */

      else {

        if (n < set_con->max_steps) {
          con->step *= 0.5;
          con->param   = con->param_old + con->step;
          assign_parameter_conwrap(con->param);

          if (order == 1) {
            for (i = 0; i < con->numUnks; i++)
              con->x[i] =  con->x_old[i] + con->step * con->x_tang[i];
          }
          else {

            con->arc_step *= 0.5;
            for (i = 0; i < con->numUnks; i++)
              con->x[i] =  con->x_old[i] + con->arc_step * con->x_tang[i];
          }
        }

        /*
         * If it was the last step, however, reset to previous solution.
         */

        else {
          con->param  = con->param_old;
          assign_parameter_conwrap(con->param);
          con->step = 0.0;
          for (i = 0; i < con->numUnks; i++)
            con->x[i] = con->x_old[i];
        }

        /*
         * Print out failure message
         */

        if (con->printproc) print_cont_step_fail(order, con);
      }
    }

    else {

      /*
       * Solver did Converge!!
       * Print out final results of a successful time step
       */

      if (order == 2) con->step = con->param - con->param_old;

      if (con->printproc) print_cont_step2(order, con);

      /*
       * If first continuation step, set to value from input file.  If
       * controlled steps, use # of Newton iters to pick next step size
       * Note:  without time step control, step size can never increase.
       */

      con->step_old = con->step;
      if (con->step_num == 0) con->step = set_con->first_step;
      else  {
	/* normal step control */
	if (set_con->step_ctrl > 0.0) {
          if (order == 2) con->arc_step *= simple_step_control(num_newt_conv,
  			         set_con->max_newton_its, set_con->step_ctrl);
          else            con->step *= simple_step_control(num_newt_conv,
  			         set_con->max_newton_its, set_con->step_ctrl);
	}
        /* for constant step runs where the step has been cut, let it
	 * increase again with step control of 0.5
	 */
        else if (order < 2 && fabs(con->step) < fabs(set_con->first_step)) {
          if (order == 2) con->arc_step *= simple_step_control(num_newt_conv,
  			         set_con->max_newton_its, 0.5);
          else            con->step *= simple_step_control(num_newt_conv,
  			         set_con->max_newton_its, 0.5);
	}
      }

      /*
       * Output information at the end of every successful time step
       * For phase transitions output second step. The Null vector for 
       * turning point calculations could be output the same way.
       */

       solution_output_conwrap(con->x, con->param, 0,
		               con->step_num, num_newt_conv);

       if (con->method == PHASE_TRANSITION_CONTINUATION)
          solution_output_conwrap(con->x_tang, con->param, 1,
		               con->step_num, num_newt_conv);

      /*
       * Check current parameter value against the maximum.
       */

      if (con->param >= set_con->max_param) {
        n = set_con->max_steps;       /* Force IO and exit  */
        if (con->printproc) {
          printf("Simulation completed continuatiion in %d steps\n",
                        con->step_num);
          printf("Final Parameter Value: %e\n\n", con->param);
        }
      }

      if (n < set_con->max_steps) {

        /*
         * Finally, its time to do some continuation, since the previous step
         * converged and it wasn't the last step.
         */

        if (con->printproc) {
          printf("\n"); print_line("~", 80);
          printf("\nCalculating initial guess for next continuation "
                        "step\n");
        }

        /*
         * Set the continuation order 0, 1 (Euler-Newton), or 2 (Arc-length)
         * Always do only zero order continuation of a turning point
	 * or any other special continuation.
         */

	switch (con->method) {
	  case FIRST_ORDER_CONTINUATION:  order = 1; break;
	  case ARC_LENGTH_CONTINUATION:  order = 2; break;
	  default: order = 0;
        }

        /*
         * Possibly adjust the con->step value for this step so it hits maximum
         * value exactly.
         */

        if ( (order < 2) &&
             (con->param + con->step >= set_con->max_param) ) {
          con->step = set_con->max_param - con->param;
        }

        /*
         * Calculate the tangent to the solution, con->x_tang, for the
	 * current step. This is trivial for 0-order continuation, requires
	 * 1 linear solve for 1st order continuation and for arc-length
	 * continuation. Use tangent to predict new solution in con->x.
         */

        switch (order) {

        case 0:
          if (con->printproc) {
            printf("\n   Doing Zeroth-order continuation --");
            printf("\n   previous solution used as initial guess\n");
          }

	  /* NO definition of con->x_tang needed for zero order.
	   * Don't set it to zero because that will mess up the
	   * turning point and phase transition tracking algorithms,
	   * which use that space for other purposes.
	   */

          /*
           * Save the old solution, before overwriting with the new solution
           */

          dcopy_(&con->numUnks, con->x, &one, con->x_old, &one);

          /* perturb guess off of singularity in tp_continuation */

          if (con->method == TURNING_POINT_CONTINUATION ||
              con->method == PITCHFORK_CONTINUATION) {
             if (con->printproc) printf("\tcon_lib: Adding random"
                                        " perturbation for tp_continuation\n");
             perturb_solution_conwrap(con->x, con->x_old,
	           	              con->scale_vec, con->numOwnedUnks);
          }

          break;

        case 1:
          if (con->printproc) {
            printf("\n   Doing First-order continuation --");
            printf("\n   calculating tangent vector by one linear "
                          "solve\n");
          }

          /* Put current solution temporarily in con->x_tang */

          dcopy_(&con->numUnks, con->x, &one, con->x_tang, &one);

          /*
           * Choose perturbation for numerical derivative of Residuals w.r.t
           * continuation parameter, and solve for the tangent vector as in
           * eq. 7.13 in John's thesis.  The continuation parameter and
           * perturbation, con->param and con->delta_param, are passed to the
	   * linear solver in the spots for the time and CJ, which aren't'
	   * needed for steady problems.
           */

          con->delta_param = con->step * 1.0e-3;

          calc_rhs_continuation(CONT_TANGENT, con->x_tang, NULL, NULL, NULL,
			        con->delta_param, con->param, NULL,
				con->numUnks, con->numOwnedUnks);

          i = linear_solver_conwrap(con->x_tang, NEW_JACOBIAN, NULL);

          /*
           * Save the old solution, before overwriting with the new solution
           */

          dcopy_(&con->numUnks, con->x, &one, con->x_old, &one);

          /*
           * Multiply the tangent vector, con->x_tang initially, by the step
           * length, and add to con->x, to obtain an initial guess at
           * the next parameter value.
           */

          for (i = 0; i < con->numUnks; i++) {
            con->x[i] -= con->x_tang[i]*con->step;
          }

          break;

        case 2:
          if (con->printproc) {
            printf("\n\tDoing Pseudo Arc-length continuation --");
            printf("\n\tCalculating tangent vector by one linear "
                                                                "solve\n");
          }

          /*
           * Start by doing the same as for 1st order continuation above.
           */

          dcopy_(&con->numUnks, con->x, &one, con->x_tang, &one);

          con->delta_param = con->step * 1.0e-3;

          calc_rhs_continuation(CONT_TANGENT, con->x_tang, NULL, NULL, NULL,
			        con->delta_param, con->param, NULL,
			        con->numUnks, con->numOwnedUnks);
          i = linear_solver_conwrap(con->x_tang, NEW_JACOBIAN, NULL);

          /* Calculate average of each variable for scaling of arc-length eq */

	  calc_scale_vec_conwrap(con->x, con->scale_vec, con->numUnks);

          /*
           * Calculate deriv of parameter w.r.t arc length, Eq.7.14a in JNS
           * thesis
           */

          con->dp_ds = 1.0 / sqrt(1.0 + 
                          scaled_dot_prod(con->x_tang, con->x_tang,
		          	          con->scale_vec, con->numOwnedUnks));

          /*
           * If this is the first step, pick con->arc_step so that this step
           * will progress the parameter by approximately con->step
           */

          if (con->step_num == 0) {

            if (con->step < 0) {
              con->dp_ds *= -1.0;
              sn_old = -1;
            }
            else
              sn_old = 1;

            con->arc_step = con->step / con->dp_ds;
          }
          else {

           /*
            * Pick sign of con->dp_ds according to eq. 7.14b in JNS thesis
            * NOTE: -1.0 factor multiplying solution terms is because
            * con->x_tang is currently the negative of the tangent vector.
            *
            * and check if a turning point was passed --
            */

            if (-1.0 * (scaled_dot_prod(con->x_tang, con->x,
					con->scale_vec, con->numOwnedUnks) -
                        scaled_dot_prod(con->x_tang, con->x_old,
				        con->scale_vec, con->numOwnedUnks)) +
                con->param - con->param_old < 0.0) {

              con->dp_ds *= -1.0;
              sn_new = -1;
            }
            else
              sn_new = 1;

            if ((con->printproc) && sn_old != sn_new)
              printf("\n\n\tA turning point was passed !!!!!!!\n");

            sn_old = sn_new;
          }

          /*
           * Save the old solution, before overwriting with the new solution
           */
          dcopy_(&con->numUnks, con->x, &one, con->x_old, &one);

          /*
           * Calculate prediction for next step from Eqs. 7.15&7.16 in JNS
           * thesis (leaving con->x_tang = u_dot).
           */

          for (i = 0; i < con->numUnks; i++) {
            con->x_tang[i]   *= -con->dp_ds;
            con->x[i] += con->x_tang[i] * con->arc_step;
          }
          con->step = con->dp_ds * con->arc_step;

          break;
        }

        /*
         * Increment the continuation parameter.  Update the
         * BCs/Properties/whatever that the continuation parameter really
         * represents.
         */

        con->param_old = con->param;
        con->param += con->step;
        assign_parameter_conwrap(con->param);

        /*
         * Increment the step counter. Print final message.
         */

        con->step_num++;

      }  /* END of:  if (n < set_con->max_steps) */

    }  /* END of else section for converged solves */

  } /* END of loop over continuation step attempts --- for (n = 0; ... --- */

  /*********************CLEAN-UP AREA*****************************************/

  /*
   * Free auxillary vectors no matter what happened
   */

  safe_free((void **) &con->x_tang);
  safe_free((void **) &con->x_old);
  safe_free((void **) &con->scale_vec);

  /*
   * Send back the overall result of the time step
   */

  return con->step_num;

} /**************** END of solve_continuation () *****************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double scaled_dot_prod(double *x, double *y,
	                      double *scale_vec, int n)
{
  int i;
  double sum=0;

  for (i=0; i<n; i++)
     sum += x[i]*y[i]*scale_vec[i]*scale_vec[i]*1.0e4;
  return gsum_double_conwrap(sum);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static double simple_step_control(int num_newt_conv, int max_Newton_steps,
  	              		  double step_ctrl)
/*
 * Function to calculate the increase in time step for the pseudo time-
 * integration option based on:
 *
 *    num_newt_conv,      the number of Newton iterations the last step
 *                        required to reach convergence, and
 *    max_Newton_steps,   the maximum number of Newton steps allowed.
 *    step_ctrl           aggressiveness of step size routine,
 *                        0.0 for constant step, 2.0 is very big
 *
 * This simplistic function will increase the time step if the last step
 * converged. It is a quadratic function in the ratio of the number of 
 * Newton steps taken compared to the maximum number of steps allowed
 * up to a maximum of 1+aggressiveness times the previous step.
 */

{

  double factor;

  factor = (max_Newton_steps - (double) num_newt_conv) /
    (max_Newton_steps - 1.0);

  return(1.0 + step_ctrl * factor * factor);

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_cont_step1(int order, struct con_struct *con)

/*
 * Print out for relevant time step information
 */

{

  char *string;

  if (con->method == TURNING_POINT_CONTINUATION)
             string = "Zero-order Turning Point Continuation";
  else if (con->method == PHASE_TRANSITION_CONTINUATION)
             string = "Zero-order Phase Transition Continuation";
  else if (con->method == PITCHFORK_CONTINUATION)
             string = "Zero-order Pitchfork Continuation";
  else if (order == 0) string = "Zero-order Continuation";
  else if (order == 1) string = "First-order Continuation";
  else if (order == 2) string = "Pseudo Arc-length Continuation";

  printf("\n"); print_line("~", 80);
  printf("\nStart of Step: %5d   Continuation Param = %9.5g from "
                "%9.5g\n", con->step_num, con->param, con->param - con->step);
  printf("\tContinuation method = %s\n", string);
  printf("\tdelta_c_p        = %8.5e ", con->step);
  printf("\n\tdelta_c_p_old    = %8.5e\n", con->step_old);

  if (order == 2) printf("\tdelta_s = %g\n", con->arc_step);
  printf("\n"); print_line("~", 80);

} /*************** END print_cont_step1 **************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_cont_step2(int order, struct con_struct *con)

/*
 * Print out for relevant continuation step information
 */

{
  printf("\n"); print_line("~", 80);
  printf("Continuation Step Number %5d was a success: Param = %10g\n",
                con->step_num, con->param);
  if (con->method == TURNING_POINT_CONTINUATION)
      printf("\tTurning Point located at:  %10g  %10g\n",
                con->param, con->tp_param);
  else if (con->method == PHASE_TRANSITION_CONTINUATION)
      printf("\tPhase Transition located at:  %10g  %10g\n",
                con->param, con->tp_param);
  else if (con->method == PITCHFORK_CONTINUATION)
      printf("\tPitchfork Bifurcation located at:  %10g  %10g\n",
                con->param, con->tp_param);
  else printf("\tOrder   = %d\n", order);

  if (order == 2) printf("\t delta_s = %g", con->arc_step);
  printf("\tdelta_c_p = %g\n", con->step);
  printf("\n"); print_line("~", 80); printf("\n");

} /************* END of print_cont_step2 () **********************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

static void print_cont_step_fail(int order, struct con_struct *con)

/*
 * Print Out descriptive information on why the current step failed
 */

{
  printf("\n"); print_line("~", 80);
  printf("\tContinuation Step Number %5d experienced a convergence "
         "failure\n", con->step_num);
  printf("\tin the non-linear or linear solver\n");
  printf("\t\tValue of parameter at failed step    = %g\n", con->param);

  if (order < 2) {
    printf("\t\tdelta_c_p of the failed step         = %g\n", 2.0* con->step);
    printf("\t\tNext value of delta_c_p              = %g\n", con->step);
  }
  else {
    printf("\t\tdelta_s of the failed step         = %g\n", 2.0*con->arc_step);
    printf("\t\tNext value of delta_s              = %g\n", con->arc_step);
  }
  printf("\n"); print_line("~", 80);

} /*************** END of print_cont_step_fail () ****************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static void print_line (char *charstr, int ntimes)

{
  int i;
  for (i = 0; i < ntimes; i++) printf("%c", *charstr);
  printf("\n");
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
