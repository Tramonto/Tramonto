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
#ifndef DFT_VEC_SPACE_H
#define DFT_VEC_SPACE_H

#include "dft_pmachine.h"

/*! \file dft_vec_space.h
\brief File defining struct that contains all information about a vector space object.

The file must be included by each Tramonto source file that needs access
to linear algebra functionality.  It contains information on the size 
and distribution of linear algebra objects such as vectors and matrices
on the parallel machine.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \def DFT_AZTEC
  \brief Selector for Aztec implementation of Tramonto's vector space interface.
*/
#define DFT_AZTEC           1

/*! \def DFT_EPETRA
  \brief Selector for Epetra's implementation of Tramonto's vector space interface.
*/
#define DFT_EPETRA          2

/*! \def DFT_LAPACK
  \brief Selector for LAPACK's implementation of Tramonto's vector space interface.
*/
#define DFT_LAPACK          3

/*! \struct DFT_VEC_SPACE_STRUCT dft_vec_space.h

\brief Contains all information about a vector space object.

*/
struct DFT_VEC_SPACE_STRUCT {

  int lin_alg_lib;  /*!< Selector for linear algebra library */
  int N_loc;        /*!< Local dimension of vector space */
  int * update;     /*!< Array of length N_loc containing GIDs owned by this processor */
  DFT_PMACHINE * machine /*!< A valid pointer to a DFT_PMACHINE object */
  void * aux_info;  /*!< Generic pointer that can be used by the implementing library for addional data */
  
};

/*! \typedef DFT_VEC_SPACE

\brief Convenient expression for struct DFT_VEC_SPACE_STRUCT.

*/
typedef struct DFT_VEC_SPACE_STRUCT DFT_VEC_SPACE;

/*! \fn void create_vec_space(int lin_alg_lib, int N_loc, int * update, DFT_PMACHINE * machine,
                              DFT_VEC_SPACE ** vec_space)

\brief Create a vector space of size N_loc on each processor with global IDs from update.

\param lin_alg_lib (In) Integer indicating the linear algebra library presently in use.
\param N_loc (In) Size of this processor's portion of the vector space.
\param update (In) Array of global IDs to use for this portion of the vector space.
\param machine (In) A valid DFT_PMACHINE pointer.
\param vec_space (Out) The address of the vector space.

\pre lin_alg_lib, N_loc, update and comm must all be properly defined.
\pre vec_space must be an uninitialized pointer.  If it is already initialized, there will be 
a memory leak.

\post The contents of the update array will be copied, so this array can be deleted without affecting
the vec_space object.  
\post The comm pointer will be copied, but the associated communicator WILL NOT be copied.  It is the
responsibility of the calling program to maintain the integrity of the associated MPI communicator.

*/

void create_vec_space(int lin_alg_lib, int N_loc, int * update, DFT_PMACHINE * machine,
                      DFT_VEC_SPACE ** vec_space);

/*! \fn void destroy_vec_space(DFT_VEC_SPACE ** vec_space)

\brief Destroy a vector space and delete all associated memory.

\param vec_space (In/Out) The address of the vector space object that will be destroyed.

\pre vec_space must have been created by a call to create_vec_space.

\post All memory associated with this object will be deleted and the pointer will be set to NULL.  
\post The associated MPI communicator will NOT be deleted.  This is the responsibility of the calling
program.

*/

void destroy_vec_space(DFT_VEC_SPACE ** vec_space);

/*! \fn void dft_random_vector(DFT_VEC_SPACE * vec_space, double *x)

\brief Fill a vector with random numbers.

\param machine (In) A pointer to a valid parallel machine object.
\param x (Out) Array of pseudo-random doubles of length compatible with the vector space.

*/

double dft_random_vector(DFT_PMACHINE * machine, double *x);

#endif /* DFT_VEC_SPACE_H */
