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
static char *cvs_util_id =
  "$Id$";
#endif
*/

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>

extern double gsum_double_conwrap(double sum);

int    N_o = -1;  /* Length of vector that is acted on (owned unknowns)  */
int    N_t = -1;  /* Length of vector that is allcoated (total unknowns) */
double N_g = -1;  /* Total number of unknowns on all procs, cast to a
                     double, which is the sum of N_o's */

/********** R O U T I N E S   I N   T H I S   F I L E  ***********************

       NAME                     TYPE                    CALL BY
----------------------------------------------------------------------
        intiialize_util_routines        void
        init_vec                        void
        vec_copy                        void
        alloc_vec                       double *
        dp                              double
        ip                              double
******************************************************************************/

/******************* PROTOTYPES FOR STATIC FUNCTIONS *************************/


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
void initialize_util_routines(int n_o, int n_t)
/*
 *  This routine sets up information that is used by the rouitnes in this file
 * and must be called before any of the routines within it are called.
 * The current input are "int n_o" the number of owned unknowns for this
 * processors, which is the length of a vectorr that is acted upton, and
 * "int n_t" which is the length that the vector must be allocated
 */
{
  N_o = n_o;
  N_t = n_t;
  N_g = gsum_double_conwrap((double)N_o);
} /* END of routine initialize_util_routines *********************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void init_vec(double *u)

/* Initialize the vector to zero.  */

{
  register int i;

  for (i = 0; i < N_t; i++) u[i] = 0.0;

} /* END of routine init_vec *************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double *alloc_vec()

/* Allocate a vector of length N_t */

{
  double *x;
  x = (double *) malloc(N_t* sizeof(double));
  init_vec(x);
  return x;

} /* END of routine alloc_vec *******************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void free_vec (double **ptr)
{
/*
 *  This version of free calls the system's free function
 *  with maximum error checking. It also doesn't call free if ptr is
 *  the NULL pointer.
 */

  if (*ptr != NULL) {

    free((void *) *ptr);

    /*
     *  Set the value of ptr to NULL, so that further references
     *  to it will be flagged.
     */

    *ptr = NULL;
  }
}  /* free_vec */


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void vec_copy(double *dx, double *dy)

/*
 * This function copies the value of the vector of type double, dx, into the
 * vector of type double, dy.  No error checking is done. If N is negative or
 * zero, dy is not changed.  A stride of 1 is assumed, unlike the linpack
 * routine called DCOPY (which is why the name is dcopy1).
 */

{
  register int i;
  for (i = 0; i < N_t; i++) dy[i] = dx[i];

} /* END of routine vec_copy *************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
double dp(double *x, double *y)
/* simple dot product */
{
  int i;
  double sum=0.0;

  for (i=0; i<N_o; i++) sum += x[i]*y[i];
  return gsum_double_conwrap(sum);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double ip(double *x, double *y)
/* inner product for pitchfork tracking. This inner product must
   have the same symmetry as the pitchfork, so a simple dot product
   only works if the mesh is symmetric. Otherwise, this product
   must be weighted by the lumped mass at each node */
{
  return dp(x,y);
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double ltransnorm(double *x, double *scale_vec)
{
  /* rescale by factor to make averages array element size one */
  return (dp(x, scale_vec) / N_g);
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
