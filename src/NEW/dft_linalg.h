/*@HEADER
// ***********************************************************************
// 
//                Tramonto: Molecular Theories Modeling Code
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// Questions? Contact Laura J.D. Frink (ljfrink@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/
#ifndef DFT_LINALG_H
#define DFT_LINALG_H

#include "dft_vec_space.h"

/*! \file dft_linalg.h
\brief Function prototypes required for Tramonto vector functionality

The file must be included by each Tramonto source file that needs access
to linear algebra functionality.  This file also specifies the functions that a 
linear algebra library must implement in order to supply vector functionality
for Tramonto.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \fn void dft_gather_global_vec(DFT_VEC_SPACE * vec_space, double *loc_vec, 
                                   int *loc_index,int N_loc, double *global_vec)

\brief Collect a distributed vector on root processor.

\param vec_space (In) A pointer to a valid vector space object.
 \param loc_vec (In) vector containing the solution vector index by the 
        local node number on the current processor.
 \param local_index (In) vector containing global node numbers of each of 
        the local processor nodes.
 \param N_loc (In) Length of loc_vec and loc_index arrays.
 \param global_vec (Output) global vector which has been gathered from all procs.

*/

void dft_gather_global_vec(DFT_VEC_SPACE * vec_space, double *loc_vec, 
                           int *loc_index,int N_loc, double *global_vec);

/*! \fn void dft_gather_global_vec_int(DFT_VEC_SPACE * vec_space, int *loc_vec, 
                                       int *loc_index,int N_loc, int *global_vec)

\brief Collect a distributed vector on root processor.

\param vec_space (In) A pointer to a valid vector space object.
 \param loc_vec (In) vector containing the solution vector index by the 
        local node number on the current processor. 
 \param local_index (In) vector containing global node numbers of each of 
        the local processor nodes.
 \param N_loc (In) Length of loc_vec and loc_index arrays.
 \param global_vec (Output) global vector which has been gathered from all procs.

*/

void dft_gather_global_vec_int(DFT_VEC_SPACE * vec_space, int *loc_vec, 
                               int *loc_index,int N_loc, int *global_vec);

/*! \fn void dft_matvec(DFT_VEC_SPACE * vec_space, void * mat, double *x, double *y)

\brief Sparse matrix-vector multiply kernel.

\param vec_space (In) A pointer to a valid vector space object.
\param A (In) A generic pointer to a sparse matrix object that will be 
       appropriately recast based on the vec_space object.
\param x (In) Array of doubles of length compatible with the vec_space.
\param y (Out) Array of doubles of length compatible with the vec_space.

*/

void dft_matvec(DFT_VEC_SPACE * vec_space, void * mat, double *x, double *y);

/*! \fn dft_gsum_double(DFT_VEC_SPACE * vec_space, double partial_sum)

\brief Global sum.

\param vec_space (In) A pointer to a valid vector space object.
\param partial_sum (In) Partial sum from local results.

\return Global sum of partial sum values.

*/

double dft_gsum_double(DFT_VEC_SPACE * vec_space, double partial_sum);

/*! \fn dft_gmax_double(DFT_VEC_SPACE * vec_space, double partial_max)

\brief Global max.

\param vec_space (In) A pointer to a valid vector space object.
\param partial_sum (In) Partial max from local results.

\return Global max of partial max values.

*/
double dft_gmax_double(DFT_VEC_SPACE * vec_space, double partial_max);


/*! \fn void dft_random_vector(DFT_VEC_SPACE * vec_space, double *x)

\brief Sparse matrix-vector multiply kernel.

\param vec_space (In) A pointer to a valid vector space object.
\param x (Out) Array of pseudo-random doubles of length compatible with the vec_space.

*/

double dft_random_vector(DFT_VEC_SPACE * vec_space, double *x);




#endif /* DFT_LINALG_H */
