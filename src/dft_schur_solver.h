/*
//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

#ifndef DFT_SCHUR_SOLVER_H
#define DFT_SCHUR_SOLVER_H

#ifdef __cplusplus
extern "C" {
#endif


/*! \file dft_schur_solver.h
\brief File defining all information needed for a Schur complement solver.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \struct DFT_SCHUR_SOLVER_STRUCT dft_schur_solver.h

\brief Contains all information about Schur complement solver.

*/
struct DFT_SCHUR_SOLVER_STRUCT {

  void * A11;  /*!< Generic pointer to 11 block matrix in Epetra_CrsMatrix form */
  void * A12;  /*!< Generic pointer to 12 block matrix in Epetra_CrsMatrix form */
  void * A21;  /*!< Generic pointer to 21 block matrix in Epetra_CrsMatrix form */
  void * A22;  /*!< Generic pointer to 22 block matrix in Epetra_CrsMatrix form */
  void * RowMap;  /*!< Generic pointer to the row map of the complete operator */
  void * ColMap;  /*!< Generic pointer to the column map of the complete operator */
  void * A;  /*!< Generic pointer to full matrix in Epetra_CrsMatrix form */
  
};

/*! \typedef DFT_SCHUR_SOLVER

\brief Convenient expression for struct DFT_SCHUR_SOLVER_STRUCT.

*/
typedef struct DFT_SCHUR_SOLVER_STRUCT DFT_SCHUR_SOLVER;

/*! \fn void dft_create_schur_solver(int * proc_config, int * external, int * bindx, 
			     double * val, int * update, int * update_index,
			     int * extern_index, int * data_org, int N_update, 
			     int Nunk_per_node,
			     DFT_SCHUR_SOLVER ** schur_solver)

\brief Create a Schur complement solver from Aztec matrix.

\param proc_config (In) Aztec array.
\param external (In) Aztec array.
\param bindx (In) Aztec array.
\param val (In) Aztec array.
\param update (In) Aztec array.
\param update_index (In) Aztec array.
\param extern_index (In) Aztec array.
\param data_org (In) Aztec array.
\param N_update (In) Aztec N_update value.
\param Nunk_per_node (In) Tramonto variable containing the number of unknowns per node on the mesh.
\param schur_solver (Out) The address of the created Schur complement solver object.

\pre All Aztec-related arguments must be properly defined, presumably by a call to AZ_transform.
\pre schur_solver must be an uninitialized pointer.  If it is already initialized, there will be 
a memory leak.

\post All input arguments are untouched except the schur_solver pointer.
\post A collection of Epetra object are created and stored in a DFT_SCHUR_SOLVER_STRUCT struct.

*/

void dft_create_schur_solver(int * proc_config, int * external, int * bindx, 
			     double * val, int * update, int * update_index,
			     int * extern_index, int * data_org, int N_update, 
			     int Nunk_per_node,
			     DFT_SCHUR_SOLVER ** schur_solver);


/*! \fn void dft_update_schur_solver(int * bindx, double * val, 
                                     DFT_SCHUR_SOLVER * schur_solver)

\brief Update the matrix values of the Schur complement operators using values in val.

\param bindx (In) Aztec bindx array.
\param val (In) Aztec val array.
\param update_index (In) Aztec update_index array.
\param schur_solver (In/Out) The address of the Schur complement solver object.

\pre bindx and val must be valid.
\pre schur_solver must be a pointer to the object created by dft_create_schur_solver, using an Aztec matrix
that is identical to the present one, with the exception of the values.  

\post The schur_solver pointer will be updated with new values from val.  bindx and val are unchanged.

*/

  void dft_update_schur_solver(int * bindx, double * val, int * update_index,
			     DFT_SCHUR_SOLVER * schur_solver);

/*! \fn void dft_apply_schur_solver(DFT_SCHUR_SOLVER * schur_solver, double * x, double * b)

\brief Solve for x given b, using the Schur complement solver.

\param vecspace (In) Vector space containing information to create operator.
*/
void dft_apply_schur_solver(DFT_SCHUR_SOLVER * schur_solver, int * update_index, double * x, double * b);


/*! \fn void dft_destroy_schur_solver(DFT_SCHUR_SOLVER ** schur_solver)

\brief Destroy a Schur complement solver object and delete all associated memory.

\param schur_solver (In/Out) The address of the linear operator object that will be destroyed.

\pre schur_solver must have been created by a call to dft_create_schur_solver.

\post All memory associated with this object will be deleted and the pointer will be set to NULL.  

*/

void dft_destroy_schur_solver(DFT_SCHUR_SOLVER ** schur_solver);
  void dft_schur_solver_comp_invA11(void * A11);

#ifdef __cplusplus
}
#endif

#endif /* DFT_SCHUR_SOLVER_H */
