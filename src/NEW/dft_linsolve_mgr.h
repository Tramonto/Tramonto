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
#ifndef DFT_LINSOLVE_MGR_H
#define DFT_LINSOLVE_MGR_H

#include "dft_vec_space.h"

/*! \file dft_linsolve_mgr.h
\brief File defining struct that contains all information to manage the linear solution process.

The file must be included by each Tramonto source file that needs access
to linear solver functionality.  It contains information on the solution strategy
and linear algebra objects that will be used to solve the linear problem.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \def DFT_AZTECSOLVE
  \brief Selector for Aztec implementation of Tramonto's linear solver interface.
*/
#define DFT_AZTECSOLVE           1

/*! \def DFT_AZTECOOSOLVE
  \brief Selector for AztecOO's implementation of Tramonto's linear solver interface.
*/
#define DFT_AZTECOOSOLVE          2

/*! \def DFT_LAPACKSOLVE
  \brief Selector for LAPACK's implementation of Tramonto's linear solver interface.
*/
#define DFT_LAPACKSOLVE          3

/*! \struct DFT_LINSOLVE_MGR_STRUCT dft_linsolve_mgr.h

\brief Contains all information about a linear solver manager object.

*/
struct DFT_LINSOLVE_MGR_STRUCT {

  int lin_solver_lib;  /*!< Selector for linear solver library */
  dft_vec_space * vec_space;  /*!< Pointer to a valid vector space object */
  void * aux_info;  /*!< Generic pointer that can be used by the implementing library for addional data */
  
};

/*! \typedef DFT_LINSOLVE_MGR

\brief Convenient expression for struct DFT_LINSOLVE_MGR_STRUCT.

*/
typedef struct DFT_LINSOLVE_STRUCT DFT_LINSOLVE_MGR;

/*! \fn void create_linsolve_mgr(int lin_solver_lib, DFT_VEC_SPACE * vec_space,
                                 DFT_LINSOLVE_MGR ** mgr)

\brief Create a linear solver manager.

\param lin_solve_lib (In) Integer indicating the linear solver library presently in use.
\param vec_space (In) A pointer to a vector space.
\param mgr (In) The address of a new linear solver manager.

\pre lin_solver_lib and vec_space must all be properly defined.
\pre mgr must be an uninitialized pointer.  If it is already initialized, there will be 
a memory leak.

\post A valid LINSOLVE_MGR object will be returned in mgr.
\post It is expected that the user will maintain the viability of the vec_space object as
long as the mgr object exists.
*/

void create_linsolve_mgr(int lin_solver_lib, DFT_VEC_SPACE * vec_space,
                      DFT_LINSOLVE_MGR ** mgr);

/*! \fn void destroy_linsolve_mgr(DFT_LINSOLVE_MGR ** mgr)

\brief Destroy a linear solver manager and delete all associated memory.

\param mgr (In/Out) The address of the linear solver manager object that will be destroyed.

\pre mgr must have been created by a call to create_linsolve_mgr.

\post All memory associated with this object will be deleted and the pointer will be set to NULL.  
\post The associated vec_space object will NOT be deleted.  This is the responsibility of the calling
program.

*/

void destroy_linsolve_mgr(DFT_LINSOLVE_MGR ** mgr);

#endif /* DFT_LINSOLVE_MGR_H */
