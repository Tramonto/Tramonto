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
extern void    init_vec(double *u);
extern void    vec_copy(double *dx, double *dy);
extern double  dp(double *x, double *y);
extern double  ip(double *x, double *y);
extern double  ltransnorm(double *x, double *y);
extern double *alloc_vec();
extern void    free_vec(double **ptr);

