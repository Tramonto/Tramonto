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

#include <stdlib.h>
#include <stdio.h>

#ifdef __STDC__
#include <stdarg.h>
#else
#include <varargs.h>
#endif

#ifdef DEBUG
    extern int Proc;
#endif

/* prototypes */
static void * smalloc(size_t n);
#ifdef __STDC__
void * array_alloc(int numdim, ...);
#else
void * array_alloc(...);
#endif

/* user callable interface to array_alloc */

void * array_alloc_1d ( size_t n1, size_t size )
{
  return array_alloc( (int) 1, (size_t) n1, (size_t) size );
}

void * array_alloc_2d ( size_t n1, size_t n2, size_t size )
{
  return array_alloc( (int) 2, (size_t) n1, (size_t) n2, (size_t) size );
}

void * array_alloc_3d ( size_t n1, size_t n2, size_t n3, size_t size )
{
  return array_alloc( (int) 3, (size_t) n1, (size_t) n2, (size_t) n3, (size_t) size );
}

void * array_alloc_4d ( size_t n1, size_t n2, size_t n3, size_t n4, size_t size )
{
  return array_alloc( (int) 4, (size_t) n1, (size_t) n2, (size_t) n3, (size_t) n4, (size_t) size );
}



/* array_alloc variable argument routine */
/******************************************************************************
 *
 *                    Dynamic Allocation of Multidimensional Arrays
 *-----------------------------------------------------------------------------
 *
 * Example Usage:
 *
 *     typedef  struct
 *       {      int     bus1;
 *              int     bus2;
 *              int     dest;
 *      }       POINT;
 *
 *      POINT    **points, corner;
 *
 *      points = (POINT **) array_alloc (2, x, y, sizeof(POINT));
 *                               ^ ^ ^
 *                               | | |
 *         number of dimensions--+ | |
 *                                 | |
 *          first dimension max----+ |
 *                                   |
 *         second dimension max------+
 *
 *         (points may be now be used as if it were declared
 *          POINT points[x][y])
 *
 *          This particular version is limited to dimensions of 3 or less.
 *
 *      corner = points[2][3]; (refer to the structure as you would any array)
 *
 *      free (points); (frees the entire structure in one fell swoop)
 *
 *****************************************************************************/
/******************************************************************************
*       The following section is a commented section containing
*       an example main code:
*******************************************************************************
*void *array_alloc();
*main()
*{
*   int ***temp;
*   int *temp2;
*   int i, j, k;
*   int il, jl, kl;
*
*   malloc_debug(2);
*   il = 2;
*   jl = 3;
*   kl = 3;
*   temp = (int ***) array_alloc(3,il,jl,kl,sizeof(int));
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) temp[i][j][k] = 1;
*      }
*   }
*
*   temp2 = (int *) malloc(10*sizeof(int));
*   for (i=0; i<10; i++) temp2[i] = 0;
*
*   for (i=0; i<il; i++) {
*      for (j=0; j<jl; j++) {
*         for (k=0; k<kl; k++) (void) printf(" %d\n", temp[i][j][k]);
*      }
*   }
*   malloc_verify();
*}
******************************************************************************/

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

#ifdef __STDC__
void *array_alloc(int numdim, ...)
#else
void  *array_alloc(va_alist)
va_dcl
#endif

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

{
  int i, j;
  struct dim {
    size_t index;  /* Number of elements in the dimension  */
    size_t total;  /* Total number of elements             */
    size_t size;   /* Size of a single element in bytes    */
    size_t off;    /* offset from beginning of array       */
  } dim[4];      /* Info about each dimension            */

#ifndef __STDC__
  int numdim;           /* Number of dimensions                 */
#endif

  size_t  total;        /* Total size of the array              */
  void  *dfield;        /* ptr to avoid lint complaints         */
  char   *field;        /* The multi-dimensional array          */
  char  **ptr;          /* Pointer offset                       */
  char   *data;         /* Data offset                          */
  va_list va;           /* Current pointer in the argument list */

#ifdef __STDC__
  va_start(va, numdim);
#else
  va_start(va);
  numdim = va_arg(va, int);
#endif

  if (numdim <= 0) {
    fprintf(stderr, "array_alloc ERROR: number of dimensions, %d, is <=0\n",
            numdim);
    exit(EXIT_FAILURE);
  }
  else if (numdim > 4) {
    fprintf(stderr, "array_alloc ERROR: number of dimensions, %d, is > 4\n",
            numdim);
    exit(EXIT_FAILURE);
  }

  dim[0].index = va_arg(va, size_t);

  if (dim[0].index <= 0) {
#ifdef DEBUG
    fprintf(stderr, "WARNING, Proc = %d: array_alloc called with first "
            "dimension <= 0, %ld\n\twill return the nil pointer\n", Proc,
            dim[0].index);
#endif
    return((double *) NULL);
  }

  dim[0].total = dim[0].index;
  dim[0].size  = sizeof(void *);
  dim[0].off   = 0;
  for (i = 1; i < numdim; i++) {
    dim[i].index = va_arg(va, size_t);
    if (dim[i].index <= 0) {
      fprintf(stderr, "WARNING: array_alloc called with dimension %d <= 0, "
              "%ld\n", i+1, dim[i].index);
      fprintf(stderr, "\twill return the nil pointer\n");
      return((double *) NULL);
    }
    dim[i].total = dim[i-1].total * dim[i].index;
    dim[i].size  = sizeof(void *);
    dim[i].off   = dim[i-1].off + dim[i-1].total * dim[i-1].size;
  }

  /* 
     last argument passed to array_alloc is provided by the sizeof 
     operator which returns a type size_t.  size_t may be an 
     unsigned int or unsigned long depending on the platform, so 
     to be safe we should read a type size_t
  */
  dim[numdim-1].size = va_arg(va, size_t);
  va_end(va);

  /* Round up the last offset value so data is properly aligned. */

  dim[numdim-1].off = dim[numdim-1].size *
    ((dim[numdim-1].off+dim[numdim-1].size-1)/dim[numdim-1].size);

  total = dim[numdim-1].off + dim[numdim-1].total * dim[numdim-1].size;

  dfield = (void *) smalloc((size_t) total);
  field  = (char *) dfield;

  for (i = 0; i < numdim - 1; i++) {
    ptr  = (char **) (field + dim[i].off);
    data = (char *) (field + dim[i+1].off);
    for (j = 0; j < dim[i].total; j++) {
      ptr[j] = data + j * dim[i+1].size * dim[i+1].index;
    }
  }

  return dfield;

} /* array_alloc */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/* Safe version of malloc.  Does not initialize memory .*/

/* Modified by Scott Hutchinson (1421) 20 January 1993 */

static void * smalloc(size_t n)

{
  void * pntr;            /* return value */
#ifndef DEBUG             /* if we're not defining DEBUG then define this to be zero */
  int Proc = 0;
#endif

  if (n < 0) {
    fprintf(stderr, "smalloc ERROR: Non-positive argument. (%d)\n", n);
    exit(EXIT_FAILURE);
  }
  else if (n == 0)
    pntr = NULL;
  else
    pntr = (void *) malloc((size_t) n);

  if (pntr == NULL && n != 0) {
    fprintf(stderr, "smalloc (Proc = %d): Out of space - number of bytes "
            "requested = %u\n", Proc, n);
    exit(EXIT_FAILURE);
  }

  return pntr;

} /* smalloc */

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void safe_free (void **ptr)
{
/*
 *  This version of free calls the system's free function
 *  with maximum error checking. It also doesn't call free if ptr is
 *  the NULL pointer.
 */
 
  if (*ptr != NULL) {
 
    free(*ptr);
 
    /*
     *  Set the value of ptr to NULL, so that further references
     *  to it will be flagged.
     */
 
    *ptr = NULL;
  }
}  /* safe_free */
/*****************************************************************************/
/*                      END of rf_allo.c                                     */
/*****************************************************************************/
