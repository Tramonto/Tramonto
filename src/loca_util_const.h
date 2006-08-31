/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER
*/

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
/*
#ifndef lint
static char *cvs_utilconsth_id =
  "$Id";
#endif
*/

/*
 *  loca_util_const.h:
 *
 *      This include file contains function prototypes for functions in
 *      rf_util.c.
 */

/*****************************************************************************/
/*     EXTERN STATEMENTS FOR GLOBAL FUNCTIONS IN rf_util.c                   */
/*****************************************************************************/

extern void    initialize_util_routines(int n_o, int n_t);
extern void    init_scale_utils(double *scale_vec);
extern void    vec_init(double *u);
extern void    vec_copy(double *dx, double *dy);
extern double  dp(double *x, double *y);
extern double  scaled_dp(double *x, double *y);
extern double  ip(double *x, double *y);
extern double  ltransnorm(double *x, double *y);
extern double *alloc_vec();
extern void    free_vec(double **ptr);
extern double  null_vector_resid(double r_val, double i_val,
                                 double *r_vec, double *i_vec, int mm_flag);
extern void    sort_by_real(int nconv, int ncv, int ldv, double *d, double *v);

