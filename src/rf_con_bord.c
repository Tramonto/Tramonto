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
#ifndef lint
static char *cvs_contbord_id =
  "$Id$";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "rf_con_const.h"

/******************************************************************************
*
*           GLOBALS DEFINED IN THIS FILE
*
*       Function                type                       Called By
*       --------------        ---------                ----------------------
*
*
******************************************************************************/

static double dp(double *x, double *y, int n);
static double ip(double *x, double *y, int n);
double ltransnorm(double *x, double *scale_vec, int n);
extern double scaled_dot_prod(double *, double *, double *, int);

#define PERTURB 1.0e-6
int AGS_option=0;

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int arc_length_bordering_alg(double *x, double *delta_x,
		             struct con_struct *con,
		             double reltol, double abstol)

/*
 * Bordering algorithm for arc-length continuation is performed.
 * Notation is from page 181 of JNS thesis. The vector "w" is
 * already calculated (eq 7-21b), and in array delta_x.
 */

{
  double *dx2, dc_p, y;
  int i, one = 1;
  double arc_length_equation = 0.0, param_update, scaled_resid;
  double al_eq_soln_contrib,al_eq_param_contrib;

 /*************************** BEGIN EXECUTION *******************************/


  /*
   * Next, "u" is calculated and stored as dx2. (eq. 7-21a) This is exactly
   * the same linear solve as a continuation predictor step.
   */

  dx2 = (double *) array_alloc(1, con->numUnks, sizeof(double));
  dcopy_(&con->numUnks, x, &one, dx2, &one);

  dc_p = PERTURB * (fabs(con->param) + fabs(con->param_old) + PERTURB);

  calc_rhs_continuation(ARC_CONT_SOL2, dx2, NULL, NULL, NULL, dc_p,
		        con->param, NULL, con->numUnks, con->numOwnedUnks);
  i = linear_solver_conwrap(dx2, OLD_JACOBIAN, NULL);

  /*
   * Calculate the arc-length equation ("g" in bordering algorithm notation)
   */

  al_eq_soln_contrib =
      scaled_dot_prod(x, con->x_tang, con->scale_vec, con->numOwnedUnks)
    - scaled_dot_prod(con->x_old,con->x_tang,con->scale_vec, con->numOwnedUnks);

  al_eq_param_contrib = con->dp_ds*(con->param - con->param_old);

  arc_length_equation = al_eq_soln_contrib + al_eq_param_contrib
                         - con->arc_step;

  /* Calculate "y", the -update to the continuation parameter (JNS eq. 7-22) */

  y = (arc_length_equation - 
       scaled_dot_prod(con->x_tang, delta_x, con->scale_vec, con->numOwnedUnks))
     / (con->dp_ds - 
       scaled_dot_prod(con->x_tang, dx2, con->scale_vec, con->numOwnedUnks));

  /*
   * y is subtracted because the minus sign in Newton's method has been
   * switched to the update term -- the same reason delta_x is subtracted
   * in update_soln
   */

  con->param -= y;
  assign_parameter_conwrap(con->param);

  /* Modify delta_x, which will be returned to the nonlinear solver (eq.7-23)*/

  for (i=0; i < con->numUnks; i++) delta_x[i] -= y * dx2[i];

  safe_free ((void **) &dx2);

  /*
   * Check whether or not the arc-length equation and continuation param
   * are converged
   */

  param_update = fabs(y) / (reltol * fabs(con->param) + abstol);
  scaled_resid = fabs(arc_length_equation)
                         / (reltol * fabs(con->arc_step) + abstol);

  if (con->printproc) {
    printf("\n\t\tFraction of arc-step due to soln, param = %g %g\n"
           ,al_eq_soln_contrib/con->arc_step,al_eq_param_contrib/con->arc_step);
    printf("\t\tarc-length scaled resid and update = %g %g\n",
                                                 scaled_resid, param_update);
  }

  if ((param_update < 1.0) && (scaled_resid < 1.0))
     return(TRUE);
  else
     return(FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int turning_point_alg(double *x, double *delta_x, struct con_struct *con,
		      double reltol, double abstol)

/*
 * Algorithm for locking on to a turning point.
 * Theory currently comes from a TEX document of Louis Romero. (AGS 1/98)
 * Lines labeled  SCALED  are additions for new scaling (section 4 of same
 * write-up). The SCALED lines can be commented out to recover unscaled version.
 */
#define SCALE_TP
{
  double *a, *b, *c, *d, *x_tmp, dt_p, y;
  double a_big, b_big, c_big, d_big, alpha;
  int i, one = 1, gnum_unks;
  double param_update, r_update, vecnorm, tmp;
  double *r_vec=con->x_tang;
  static int first=TRUE;

/* This flag can be 0 or 1 and effects scaling */
/* Seems to work better as 0 for TP calcs      */
AGS_option=1;

 /*************************** BEGIN EXECUTION *******************************/

  /* Allocate arrays for turning point bordering algorithm */

  a = (double *) array_alloc(1, con->numUnks, sizeof(double));
  b = (double *) array_alloc(1, con->numUnks, sizeof(double));
  c = (double *) array_alloc(1, con->numUnks, sizeof(double));
  d = (double *) array_alloc(1, con->numUnks, sizeof(double));
  x_tmp = (double *) array_alloc(1, con->numUnks, sizeof(double));

  /* construct "a" vector from delta_x */

  for (i=0; i < con->numUnks; i++) a[i] = -delta_x[i];

  /*
   * Next, "b" is calculated. This is exactly
   * the same linear solve as a continuation predictor step except
   * the tp_param is perturbed, not the continuation param.
   */

  dcopy_(&con->numUnks, x, &one, b, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL2, b, NULL, NULL, NULL, dt_p,
		        con->tp_param, NULL, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(b, OLD_JACOBIAN, NULL);

  /* first time through this routine, initialize r_vec with b */

  if (first) {

    /* calculate variable averages, to be used in ltransnorm calls below */

    calc_scale_vec_conwrap(x, con->scale_vec, con->numUnks);

    vecnorm = ltransnorm(b, con->scale_vec, con->numOwnedUnks);

    for (i=0; i < con->numUnks; i++) r_vec[i] = b[i] / vecnorm;

    first = FALSE;
  }
  else {

    /* This section is optional, changing scalng vector betweem Newton iters*/

    /* calculate variable averages, to be used in ltransnorm calls below */
    calc_scale_vec_conwrap(x, con->scale_vec, con->numUnks);

    vecnorm = ltransnorm(r_vec, con->scale_vec, con->numOwnedUnks);
    if (con->printproc)
       printf("\n\t\tRescaling r_vec by %g to make it length 1\n", vecnorm);

    for (i=0; i < con->numUnks; i++) r_vec[i] /= vecnorm;
  }

  /* Rescale a and b vectors as in Louie's write-up, section 4 */

#ifdef SCALE_TP
  b_big = 1.0 / ltransnorm(b, con->scale_vec, con->numOwnedUnks); /*SCALED*/
  a_big = -ltransnorm(a, con->scale_vec, con->numOwnedUnks) * b_big; /*SCALED*/
  for (i=0; i < con->numUnks; i++) a[i] += a_big*b[i]; /*SCALED*/
  for (i=0; i < con->numUnks; i++) b[i] *= b_big; /*SCALED*/
#endif

  /* Next, "c" is calculated as a function of a and r_vec. */

  dcopy_(&con->numUnks, x, &one, c, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL3, c, a, con->scale_vec, x_tmp, dt_p,
		        con->tp_param, r_vec, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(c, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  /* Next, "d" is calculated as a function of b and r_vec. */

  dcopy_(&con->numUnks, x, &one, d, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL4, d, b, con->scale_vec, x_tmp, dt_p,
		        con->tp_param, r_vec, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(d, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  /*
   * Calculate the updates to tp_param (stored in dt_p),
   * r_vec (stored in c), and x (stored in -delta_x).
   */

  /* Rescale c and d vectors as in Louie's write-up, section 4 */

#ifdef SCALE_TP
  d_big = 1.0 / ltransnorm(d, con->scale_vec, con->numOwnedUnks); /*SCALED*/
  c_big = -ltransnorm(c, con->scale_vec, con->numOwnedUnks) * d_big; /*SCALED*/
  for (i=0; i < con->numUnks; i++) c[i] += c_big*d[i]; /*SCALED*/
  for (i=0; i < con->numUnks; i++) d[i] *= d_big; /*SCALED*/
#endif

  dt_p = ( (1.0 - 
            AGS_option*ltransnorm(r_vec, con->scale_vec, con->numOwnedUnks))
		    - ltransnorm(c, con->scale_vec, con->numOwnedUnks) )
         / ltransnorm(d, con->scale_vec, con->numOwnedUnks);

  for (i=0; i < con->numUnks; i++) 
    c[i] += dt_p * d[i] + (AGS_option-1)*r_vec[i];

  /* dt_p change meaning from Beta to alpha here */

#ifdef SCALE_TP
  dt_p = dt_p * d_big + c_big; /*SCALED*/
#endif

  for (i=0; i < con->numUnks; i++) delta_x[i] = -a[i];

  for (i=0; i < con->numUnks; i++) delta_x[i] -= dt_p * b[i];

  /*dt_p change meaning from alpha to dt_p here */

#ifdef SCALE_TP
  dt_p = a_big + dt_p * b_big; /*SCALED*/
#endif

  /* If delta_r is not nearly zero, truncate x and r steps for safety */

  if (fabs(ltransnorm(c, con->scale_vec, con->numOwnedUnks)) > 0.00001)  {
    vecnorm = 0.00001/fabs(ltransnorm(c, con->scale_vec, con->numOwnedUnks));
    if (con->printproc)
      printf("Roundoff causing nonzero delta_r, rescaling updates by %g\n",
              vecnorm);
    for (i=0; i < con->numUnks; i++) {
      delta_x[i] *= vecnorm;
      c[i] *= vecnorm;
    }
  }

  /*
   * Update the tp_param and the null vector r_vec
   */

  con->tp_param += dt_p;
  assign_tp_parameter_conwrap(con->tp_param);

  for (i=0; i < con->numUnks; i++) printf("%d.  %g  %g\n",i,r_vec[i],c[i]);
  for (i=0; i < con->numUnks; i++) r_vec[i] += c[i];

  /*
   * Check whether or not the arc-length equation and continuation param
   * are converged
   */

  param_update = fabs(dt_p) / (reltol * fabs(con->tp_param) + abstol);

  /* get update norm of r_vec, to check for convergence */

  r_update = 0.0;
  for (i=0; i < con->numUnks; i++) {
    tmp       = c[i]/(fabs(r_vec[i])*reltol + 100.0*abstol);
    r_update += tmp*tmp;
  }
  gnum_unks = gsum_double_conwrap(con->numOwnedUnks);
  r_update = sqrt( gsum_double_conwrap(r_update)/gnum_unks );

  if (con->printproc) {
    printf("\n\ttp_parameter (%g) unscaled and scaled update = %g %g\n"
                        ,con->tp_param, dt_p, param_update);
    printf("\n\tNull vector scaled update = %g\n", r_update);
  }

  safe_free ((void **) &a);
  safe_free ((void **) &b);
  safe_free ((void **) &c);
  safe_free ((void **) &d);
  safe_free ((void **) &x_tmp);

  /* return convergence status of the parameter and null vector */

  if ((param_update < 1.0) && (r_update < 10.0))                       
     return(TRUE);
  else
     return(FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int pitchfork_alg(double *x, double *delta_x, struct con_struct *con,
	          double reltol, double abstol)

/*
 * Algorithm for locking on to a Pitchfork bifurcation
 * This assumes that con->x_tang contains an asymmetrtic vector
 * the first time through, and an initial guess for the 
 * null vector every time (with both meanings the first time).
 */

{
  double *a, *b, *c, *d, *e, *f, *x_tmp, dt_p, y;
  double a_big, b_big, c_big, d_big, alpha;
  int i, one = 1, gnum_unks, n=con->numOwnedUnks;
  double param_update, r_update, vecnorm, tmp;
  double *r_vec=con->x_tang;
  double eps_update, d_eps, ipa, ipb, ipc, ipx, ipdx, ltc, ltd, lte, ltf, ltn;
  double t1, t2, t3, t4;
  static int first=TRUE;
  static double *psi, eps;

/* This flag can be 0 or 1 and effects scaling */
/* Seems to work better as 1 for PF calcs      */
AGS_option=1;

 /*************************** BEGIN EXECUTION *******************************/

  /* first time through this routine, set psi asymmetric vector to con-tang
     and set epsilon "asymmetry variable" to 0.0 */

  if (first) {
if (con->printproc) printf("\n\tIn pitchfork alg, AGS_option=%d\n",AGS_option);
    psi = (double *) array_alloc(1, con->numUnks, sizeof(double));
    eps = 0.0;
    tmp = sqrt(ip(con->x_tang, con->x_tang, n));
    if (tmp==0.0) {
      if (con->printproc) {
        printf("ERROR in pitchfork alg: Asymmetric vector must be supllied\n");
      }
      exit(-1);
    }
    for (i=0; i < con->numUnks; i++) psi[i] = con->x_tang[i]/tmp;
    first = FALSE;

    t1 = ip(x, psi, n);
    t2 = ip(x, x, n);
    t3 = ip(psi, psi, n);
    t4 = t1 / sqrt(t2*t3);
    if (con->printproc) {
      printf("\tPitchfork Alg: On input <x,psi> = %g; Should be 0.0.\n",t4);
    }
  }

  /* calculate variable averages, to be used in ltransnorm calls below */

  calc_scale_vec_conwrap(x, con->scale_vec, con->numUnks);

  /* Make sure the null vector r_vec is length 1 to begin with */

  vecnorm = ltransnorm(r_vec, con->scale_vec, n);
  if (con->printproc)
     printf("\n\t\tRescaling r_vec by %g to make it length 1\n", vecnorm);

  for (i=0; i < con->numUnks; i++) r_vec[i] /= vecnorm;

  /* Allocate arrays for turning point bordering algorithm */

  a = (double *) array_alloc(1, con->numUnks, sizeof(double));
  b = (double *) array_alloc(1, con->numUnks, sizeof(double));
  c = (double *) array_alloc(1, con->numUnks, sizeof(double));
  d = (double *) array_alloc(1, con->numUnks, sizeof(double));
  e = (double *) array_alloc(1, con->numUnks, sizeof(double));
  f = (double *) array_alloc(1, con->numUnks, sizeof(double));
  x_tmp = (double *) array_alloc(1, con->numUnks, sizeof(double));

  /* Begin real stuff */
  /* construct "a" vector from delta_x */

  for (i=0; i < con->numUnks; i++) a[i] = -delta_x[i];

  /*
   * Next, "b" is calculated. This is exactly
   * the same linear solve as a continuation predictor step except
   * the tp_param is perturbed, not the continuation param.
   */

  dcopy_(&con->numUnks, x, &one, b, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL2, b, NULL, NULL, NULL, dt_p,
		        con->tp_param, NULL, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(b, OLD_JACOBIAN, NULL);

  /* Next, "c" is calculated using just psi as rhs */

  for (i=0; i < con->numUnks; i++) c[i] = -psi[i];

  i = linear_solver_conwrap(c, OLD_JACOBIAN, NULL);

  /* Next, "d" is calculated as a function of a and r_vec. */

  dcopy_(&con->numUnks, x, &one, d, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL3, d, a, con->scale_vec, x_tmp, dt_p,
		        con->tp_param, r_vec, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(d, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  if (AGS_option==1) for (i=0; i < con->numUnks; i++) d[i] += r_vec[i];

  /* Next, "e" is calculated as a function of b and r_vec. */

  dcopy_(&con->numUnks, x, &one, e, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL4, e, b, con->scale_vec, x_tmp, dt_p,
		        con->tp_param, r_vec, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(e, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  /* Next, "f" is calculated as a function of c and r_vec. */

  dcopy_(&con->numUnks, x, &one, f, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL3, f, c, con->scale_vec, x_tmp, dt_p,
		        con->tp_param, r_vec, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(f, SAME_BUT_UNSCALED_JACOBIAN, x_tmp);

  if (AGS_option==1) for (i=0; i < con->numUnks; i++) f[i] += r_vec[i];

  /*
   * Calculate the updates to tp_param (stored in dt_p),
   * r_vec (stored in d), and x (stored in -delta_x).
   */

  ipx = ip(x, psi, n);
  ipa = ip(a, psi, n);
  ipb = ip(b, psi, n);
  ipc = ip(c, psi, n);
  ltd = ltransnorm(d, con->scale_vec, n);
  lte = ltransnorm(e, con->scale_vec, n);
  ltf = ltransnorm(f, con->scale_vec, n);
  ltn = ltransnorm(r_vec, con->scale_vec, n);

  d_eps = -eps + ( (ipx+ipa)*lte + ipb*(1.0 - ltd) ) / (ipb*ltf - ipc*lte);

  dt_p = ( 1.0 - ltd - ltf*(d_eps + eps)) / lte;

/*
  if (con->printproc)
     printf("\tPf: ipx=%g, ipa=%g, ipb=%g, ipc=%g,\n",ipx,ipa,ipb,ipc);
  if (con->printproc)
     printf("\tPf: ltd=%g, lte=%g, ltf=%g, ltn=%g\n",ltd,lte,ltf, ltn);

  if (con->printproc)
     printf("\tPitchfork info: dt_p=%g, d_eps=%g,\n", dt_p, d_eps);
*/

  /* Negative of update vector here */
  for (i=0; i < con->numUnks; i++)
     delta_x[i] = -a[i] - c[i]*(eps + d_eps) - b[i]*dt_p;

  /* use c space for delta r_vec */

  for (i=0; i < con->numUnks; i++)
    c[i] = -r_vec[i] + d[i] + e[i]*dt_p + f[i]*(eps + d_eps);

  ipdx = ip(delta_x, psi, n);
  ltc = ltransnorm(c, con->scale_vec, n);
  /*
  if (con->printproc)
     printf("ip(delta_x,psi) = %g, lt delta_n = %g, l\n", -ipdx, ltc);
  */

  /*
   * Update the tp_param, eps,  and the null vector r_vec
   */

/*
if (con->printproc) printf("writing x, n, psi, delta_x, delta_n\n");
  solution_output_conwrap(x,       1, 0, 5, 0);
  solution_output_conwrap(r_vec,   1, 0, 6, 0);
  solution_output_conwrap(psi,     1, 0, 7, 0);
  solution_output_conwrap(delta_x, 1, 0, 8, 0);
  solution_output_conwrap(c,       1, 0, 9, 0);
*/

  eps += d_eps;

  con->tp_param += dt_p;
  assign_tp_parameter_conwrap(con->tp_param);

  for (i=0; i < con->numUnks; i++) r_vec[i] += c[i];

  /*
   * Check whether or not the continuation param, null vector, and eps
   * are converged
   */

  param_update = fabs(dt_p)  / (reltol * fabs(con->tp_param) + abstol);
  eps_update   = fabs(d_eps) / (reltol * fabs(eps) + abstol);

  /* get update norm of r_vec, to check for convergence */

  r_update = 0.0;
  for (i=0; i < con->numUnks; i++) {
    tmp       = c[i]/(fabs(r_vec[i])*reltol + abstol);
    r_update += tmp*tmp;
  }
  gnum_unks = gsum_double_conwrap(con->numOwnedUnks);
  r_update = sqrt(gsum_double_conwrap(r_update)/gnum_unks);

  if (con->printproc) {
    printf("\n\ttp_parameter (%g) unscaled and scaled update = %g %g\n"
                        ,con->tp_param, dt_p, param_update);
    printf("\n\teps (%g) unscaled and scaled update = %g %g\n"
                        ,eps, d_eps, eps_update);
    printf("\n\tNull vector scaled update = %g\n", r_update);
  }

  safe_free ((void **) &a);
  safe_free ((void **) &b);
  safe_free ((void **) &c);
  safe_free ((void **) &d);
  safe_free ((void **) &e);
  safe_free ((void **) &f);
  safe_free ((void **) &x_tmp);

  /* return convergence status of the parameter, eps and null vector */

  if (param_update < 1.0 && eps_update < 1.0 && r_update < 10.0)
     return(TRUE);
  else
     return(FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int phase_transition_alg(double *x, double *delta_x,
		         struct con_struct *con, double reltol, double abstol)

/*
 * Algorithm for locking on to a phase transition.
 * This involves finding two solutions that exist at the
 * same parameters and at the same energy.
 */

{
  double *a, *b, *c, *d, *x2, g, dg_dtp, dt_p;
  double *x_tmp, *x2_tmp, dgac, dgbd, eps;
  double tmp, x2_update, param_update;
  int i, one = 1, n, gnum_unks;
  extern double free_energy_diff_conwrap(double *, double *);
  extern int update_solution(int,double *, double *);
  int converged=FALSE;

 /*************************** BEGIN EXECUTION *******************************/

  n = con->numOwnedUnks;

  x2 = con->x_tang;
  a = (double *) array_alloc(1, con->numUnks, sizeof(double));
  b = (double *) array_alloc(1, con->numUnks, sizeof(double));
  c = (double *) array_alloc(1, con->numUnks, sizeof(double));
  d = (double *) array_alloc(1, con->numUnks, sizeof(double));
  x_tmp  = (double *) array_alloc(1, con->numUnks, sizeof(double));
  x2_tmp = (double *) array_alloc(1, con->numUnks, sizeof(double));

  for (i=0; i<con->numUnks; i++) {
     a[i] = -delta_x[i];
     b[i] = 0.0;
     c[i] = 0.0;
     d[i] = 0.0;
  }

  /* Construct b vector using tangent calculation */

  dcopy_(&con->numUnks, x, &one, b, &one);

  dt_p = PERTURB * (con->tp_param + PERTURB);

  calc_rhs_continuation(TP_CONT_SOL2, b, NULL, NULL, NULL, dt_p,
		        con->tp_param, NULL, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(b, OLD_JACOBIAN, NULL);

  /* Now do Newton iteration on second vector for c */

  matrix_residual_fill_conwrap(x2, c, RHS_MATRIX);

  i = linear_solver_conwrap(c, NEW_JACOBIAN, NULL);

  for (i=0; i < con->numUnks; i++) c[i] *= -1.0;

  /* Construct d vector using tangent calculation */

  dcopy_(&con->numUnks, x2, &one, d, &one);

  calc_rhs_continuation(TP_CONT_SOL2, d, NULL, NULL, NULL, dt_p,
		        con->tp_param, NULL, con->numUnks, con->numOwnedUnks);

  i = linear_solver_conwrap(d, OLD_JACOBIAN, NULL);

  /* Now start bordering alg for equal-energy constraint
   * g is energy difference to be driven to zero,
   * dg_dtp is the derivative of the energy w.r.t. the parameter
   * dgac is the directional derivative in of g the a:c direction
   * dgbd is the directional derivative of g in the b:d direction
   */

  g = free_energy_diff_conwrap(x, x2);

  assign_tp_parameter_conwrap(con->tp_param + dt_p);

  dg_dtp = (free_energy_diff_conwrap(x, x2) - g ) / dt_p;

  assign_tp_parameter_conwrap(con->tp_param);

  eps = PERTURB * sqrt( dp(x,x,n) / (dp(a,a,n) + 1.0e-6)
                      + dp(x2,x2,n) / (dp(c,c,n) + 1.0e-6));

  for (i=0; i < con->numUnks; i++) {
    x_tmp[i]  = x[i]  + eps*a[i];
    x2_tmp[i] = x2[i] + eps*c[i];
  }
  dgac = (free_energy_diff_conwrap(x_tmp, x2_tmp) - g ) / eps;

  eps = PERTURB * sqrt(dp(x,x,n) / (dp(b,b,n) + 1.0e-6)
                     + dp(x2,x2,n) / (dp(d,d,n) + 1.0e-6));

  for (i=0; i < con->numUnks; i++) {
    x_tmp[i]  = x[i]  + eps*b[i];
    x2_tmp[i] = x2[i] + eps*d[i];
  }
  dgbd = (free_energy_diff_conwrap(x_tmp, x2_tmp) - g ) / eps;

  dt_p =  -(g  + dgac) / (dgbd + dg_dtp);

  /* update continuation parameter */

  con->tp_param += dt_p;
  assign_tp_parameter_conwrap(con->tp_param);

  /* update x2, checking for convergence */

  for (i=0; i<con->numUnks; i++) delta_x[i] = c[i] + dt_p*d[i];

  /* get update norm to check for convergence */

  x2_update = 0.0;
  for (i=0; i < con->numUnks; i++) {
    tmp      += delta_x[i]/(fabs(x2[i])*reltol + abstol);
    x2_update += tmp*tmp;
  }
  gnum_unks = gsum_double_conwrap(con->numOwnedUnks);
  x2_update  = sqrt( gsum_double_conwrap(x2_update/gnum_unks) );

  param_update = fabs(dt_p) / (reltol * fabs(con->tp_param) + abstol);

  /* update x2 and modify delta_x for update in usual Newton routine */

  /*for (i=0; i<con->numUnks; i++) x2[i]     += delta_x[i];*/
  converged=update_solution(2,x2,delta_x);
  for (i=0; i<con->numUnks; i++) delta_x[i] = -a[i] - dt_p*b[i];

  safe_free((void *) &a);
  safe_free((void *) &b);
  safe_free((void *) &c);
  safe_free((void *) &d);
  safe_free((void *) &x_tmp);
  safe_free((void *) &x2_tmp);

  if (con->printproc) {
    printf("\tPhase transition algorithm scaled parameter update  : %g\n", param_update);
    printf("\tPhase transition algorithm scaled solution #2 update: %g\n", x2_update);
  }

  /* x2_update and param_update must be < 1 for convergence */

  if (converged /*x2_update < 1.0*/ && param_update < 1.0)
     return(TRUE);
  else
     return(FALSE);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static double dp(double *x, double *y, int n)
/* simple dot product */
{
  int i;
  double sum=0.0;

  for (i=0; i<n; i++) sum += x[i]*y[i];
  return gsum_double_conwrap(sum);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
static double ip(double *x, double *y, int n)
/* inner product for pitchfork tracking. This inner product must 
   have the same symmetry as the pitchfork, so a simple dot product
   only works if the mesh is symmetric. Otherwise, this product
   must be weighted by the lumped mass at each node */
{
  int i;
  double sum=0.0;

  for (i=0; i<n; i++) sum += x[i]*y[i];
  return gsum_double_conwrap(sum);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double ltransnorm(double *x, double *scale_vec, int n)
{
  int i;
  double sum=0;

  for (i=0; i<n; i++)  sum += x[i]*scale_vec[i];

  /* rescale by factor to make averages array element size one */
  return (gsum_double_conwrap(sum)/ gsum_double_conwrap((double) n));
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void calc_rhs_continuation(int rhs_type, double *x, double *ab_vec,
		           double *scale_vec, double *x_tmp, double dc_p,
			   double param, double *r_vec,
			   int numUnks, int numOwnedUnks)

/* routine to pre-calculate the non-standard right-hand-sides for solves */
/* of continuation runs. This routine returns the matrix fill time.      */
/* The matrix at the current dolution is also calculated                 */

{
  double dc_p1, *resid_delta, *resid_vector;
  int i;

  /* Allocate and initialize resid_delta */

  resid_delta = (double *) array_alloc(1, numUnks, sizeof(double));
  resid_vector= (double *) array_alloc(1, numUnks, sizeof(double));

  for (i=0; i< numUnks; i++) {
    resid_delta[i] = 0.0;
    resid_vector[i]= 0.0;
  }

  /*
   * For the first 2 options, resid_delta is the right hand side at the
   * perturbed value of the continuation parameter
   */

  if (rhs_type == CONT_TANGENT) {

    assign_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY);

    assign_parameter_conwrap(param);

    matrix_residual_fill_conwrap(x, resid_vector, RHS_MATRIX);

    for (i=0; i< numUnks; i++)
      resid_vector[i] = (resid_delta[i] - resid_vector[i])/ dc_p;
  }

  else if (rhs_type == ARC_CONT_SOL2) {

    assign_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY);

    assign_parameter_conwrap(param);

    matrix_residual_fill_conwrap(x, resid_vector, RHS_ONLY);

    for (i=0; i< numUnks; i++)
      resid_vector[i] = (resid_delta[i] - resid_vector[i])/ dc_p;
  }

  /*
   * For TP_CONT_SOL2, resid_delta is the right hand side at the
   * perturbed value of the turning point parameter
   */

  else if (rhs_type == TP_CONT_SOL2 ) {

    assign_tp_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, RHS_ONLY);

    assign_tp_parameter_conwrap(param);

    matrix_residual_fill_conwrap(x, resid_vector, RHS_ONLY);

    for (i=0; i< numUnks; i++)
      resid_vector[i] = -(resid_delta[i] - resid_vector[i])/ dc_p;
  }

  /*
   * For TP_CONT_SOL3, resid_delta is the matrix at perturbed value of
   * solution in the direction of ab_vec multiplied by the vector r_vec
   */

  else if (rhs_type == TP_CONT_SOL3 ) {

    dc_p1 = PERTURB * (PERTURB +
		 sqrt(dp(x,x,numOwnedUnks)/(dp(ab_vec,ab_vec,numOwnedUnks)+PERTURB)));

printf("TP3: x=%g, ab=%g, del=%g\n\n",dp(x,x,numOwnedUnks),dp(ab_vec,ab_vec,numOwnedUnks),dc_p1);
    for (i=0; i< numUnks; i++) x_tmp[i] = x[i] + dc_p1 * ab_vec[i];

    matrix_residual_fill_conwrap(x_tmp, resid_delta, RHS_MATRIX);

    matvec_mult_conwrap(r_vec, resid_delta);

    matrix_residual_fill_conwrap(x, resid_vector, MATRIX_ONLY);

    matvec_mult_conwrap(r_vec, resid_vector);

/*
printf("In SOL3: n.j.n = %g  %g\n", dp(r_vec, resid_vector,numOwnedUnks),
dp(r_vec, resid_vector,numOwnedUnks)/dp(r_vec, r_vec, numOwnedUnks));
*/

    for (i=0; i< numUnks; i++)
      resid_vector[i] = - AGS_option*resid_vector[i]
                        -(resid_delta[i] - resid_vector[i])/ dc_p1;
  }

  /*
   * For TP_CONT_SOL4, resid_delta is the matrix at perturbed value of
   * solution in the direction of ab_vec multiplied by the vector r_vec
   * plus the matrix perturbed in the direction of the tp_parameter
   * multiplied by the vector r_vec.
   */

  else if (rhs_type == TP_CONT_SOL4 ) {

    dc_p1 = PERTURB * (PERTURB +
		 sqrt(dp(x,x,numOwnedUnks)/dp(ab_vec,ab_vec,numOwnedUnks)));

printf("TP4: x=%g, ab=%g, del=%g\n\n",dp(x,x,numOwnedUnks),dp(ab_vec,ab_vec,numOwnedUnks),dc_p1);

    for (i=0; i< numUnks; i++) x_tmp[i] = x[i] + dc_p1 * ab_vec[i];

    matrix_residual_fill_conwrap(x_tmp, resid_delta, MATRIX_ONLY);

    matvec_mult_conwrap(r_vec, x_tmp);

    assign_tp_parameter_conwrap(param + dc_p);

    matrix_residual_fill_conwrap(x, resid_delta, MATRIX_ONLY);

    assign_tp_parameter_conwrap(param);

    matvec_mult_conwrap(r_vec, resid_delta);

    /*
     * Two numerical differences with different perturbations are summed
     * together
     */

    for (i=0; i< numUnks; i++) resid_delta[i] =
                      resid_delta[i] / dc_p  +  x_tmp[i] / dc_p1;

    /* also sum together the perturbations */

    dc_p = 1.0 / (1.0 / dc_p1 + 1.0 / dc_p);

    matrix_residual_fill_conwrap(x, resid_vector, MATRIX_ONLY);

    matvec_mult_conwrap(r_vec, resid_vector);

    for (i=0; i< numUnks; i++)
      resid_vector[i] = -resid_delta[i] + resid_vector[i]/ dc_p;

  }

  for (i=0; i< numUnks; i++) x[i] = resid_vector[i];

  safe_free((void **) &resid_delta);
  safe_free((void **) &resid_vector);

  return;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
int continuation_hook(double *x, double *delta_x, void *con_void,
                      double reltol, double abstol)
{
  /* return TRUE for non-continuation probelms, as well as zero and
   * first order continuation runs
   */

   int converged=TRUE;
   struct con_struct *con = (struct con_struct *) con_void;

  /*
   * If arc_length continuation is being done, then branch here
   * to do the rest of the bordering algorithm. This involves a
   * second matrix fill and linear solve. (See JNS thesis
   * eqns 7.18-7.19 for equations, 7.21-7.23 for method.
   */

   if (con->method == ARC_LENGTH_CONTINUATION && con->step_num > 0)
    converged = arc_length_bordering_alg(x, delta_x, con, reltol, abstol);

  /*
   * If turning point calculations are being done, then branch here
   * to do the rest of the algorithm. This involves 6 more matrix fills
   * and 3 more linear solves in the current implementation.
   */

   else if (con->method == TURNING_POINT_CONTINUATION)
     converged = turning_point_alg(x, delta_x, con, reltol, abstol);

  /*
   * If pitchfork bifurcation calculations are being done, then branch
   * here to do the rest of the algorithm.
   */

   else if (con->method == PITCHFORK_CONTINUATION)
     converged = pitchfork_alg(x, delta_x, con, reltol, abstol);

  /*
   * If phase_transition tracking is being done, then branch here
   * to do the rest of the algorithm. This involves
   * 3 more linear solves in the current implementation.
   */

   else if (con->method == PHASE_TRANSITION_CONTINUATION)
     converged = phase_transition_alg(x, delta_x, con, reltol, abstol);

   /* No boredering alg needed for zero and first order continuations */

   else if (con->method == ZERO_ORDER_CONTINUATION ||
            con->method == FIRST_ORDER_CONTINUATION ||
	    (con->method == ARC_LENGTH_CONTINUATION && con->step_num == 0)) {
     converged = TRUE;
   }

   /* perform error check */

   else {
     if (con->printproc)
       printf("\nERROR continuation_hook: Unknown method %d\n",
		                                    con->method);
     exit(-1);
   }

   return converged;

}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
