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
#ifndef DFT_LINOP_H
#define DFT_LINOP_H

#include "dft_vec_space.h"

/*! \file dft_linop.h
\brief File defining struct that contains all information about a linear operator object.

The file must be included by each Tramonto source file that needs access
to linear operators.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \def DFT_AZTEC_LINOP
  \brief Selector for Aztec implementation of Tramonto's linear operator.
*/
#define DFT_AZTEC_LINOP           1

/*! \def DFT_EPETRA_LINOP
  \brief Selector for Epetra's implementation of Tramonto's linear operator.
*/
#define DFT_EPETRA_LINOP          2

/*! \def DFT_LAPACK_LINOP
  \brief Selector for LAPACK's implementation of Tramonto's linear operator.
*/
#define DFT_LAPACK_LINOP          3

/*! \struct DFT_LINOP_STRUCT dft_linop.h

\brief Contains all information about a linear operator object.

*/
struct DFT_LINOP_STRUCT {

  dft_vec_space * vec_space;    /*!< Selector for linear algebra library */
  int first_time_graph_insert=TRUE;   /*!< Initialized to TRUE, set to FALSE after the first call to insert_global_graph_row */
  int * internal_index; /*!< After fill_complete is called, contains the operator ordering for GIDs owned by this processor */
  int * external_index; /*!< After fill_complete is called, contains the operator ordering for GIDs not owned by this processor but present as column indices on this processor */
  void * aux_info;  /*!< Generic pointer that can be used by the implementing library for addional data */
  
};

/*! \typedef DFT_LINOP

\brief Convenient expression for struct DFT_LINOP_STRUCT.

*/
typedef struct DFT_LINOP_STRUCT DFT_LINOP;

/*! \fn void dft_create_linop(DFT_VEC_SPACE * vec_space, DFT_LINOP ** op)

\brief Create a vector space of size N_loc on each processor with global IDs from update.

\param vec_space (In) Vector space containing information to create operator.
\param op (Out) The address of the linear operator.

\pre vec_space must be properly defined.
\pre op must be an uninitialized pointer.  If it is already initialized, there will be 
a memory leak.

\post The vec_space pointer will be copied, but the associated vector space data WILL NOT be copied.  It is the
responsibility of the calling program to maintain the integrity of the associated vector space object.

*/

void dft_create_linop(DFT_VEC_SPACE * vec_space, DFT_LINOP ** op);

/*! \fn void dft_destroy_linop(DFT_LINOP ** op)

\brief Destroy a linear operator and delete all associated memory.

\param op (In/Out) The address of the linear operator object that will be destroyed.

\pre op must have been created by a call to dft_create_linop.

\post All memory associated with this object will be deleted and the pointer will be set to NULL.  
\post The associated vec_space object will NOT be deleted.  This is the responsibility of the calling
program.

*/

void dft_destroy_linop(DFT_LINOP ** op);

/*! \fn int dft_insert_global_graph_row(DFT_LINOP * op, int i_box, int i_mat, int * col_indices)

\brief Insert box column indices for i_box row into i_mat row of matrix graph.

This function passes column indices for the graph of the matrix to the DFT_LINOP object.  The column 
indices are in box coordinate space.

\param op (In) The address of the linear operator object that will be updated.

\param i_box (In) Index of row in "box" coordinate space to which these column indices belong.

\param i_mat (In) Index of row in local "matrix" coordinate space to which these column indices belong.

\param col_indices (In/Out) Dense integer array filled with TRUE/FALSE values.  col_indices[j_box] is TRUE
       if the i_box equation depends on the j_box unknown.  On exit, all values should be set to FALSE.

\pre op must have been created by a call to create_linop.

\post The graph associated with this DFT_LINOP object will be updated with the passed-in column entries.
\post This row of the matrix may not receive additional column indices and the number of column indices is fixed.

\return Returns 0 if no errors detected.

*/
int dft_insert_global_graph_row(DFT_LINOP * op, int i_box, int i_mat, int * row_indices);

/*! \fn int dft_fill_complete(DFT_LINOP * op)

\brief Insert box column indices for i_box row into i_mat row of matrix graph.

This function passes column indices for the graph of the matrix to the DFT_LINOP object.  The column 
indices are in box coordinate space.

\param op (In) The address of the linear operator object that will be updated.

\param i_box (In) Index of row in "box" coordinate space to which these column indices belong.

\param i_mat (In) Index of row in local "matrix" coordinate space to which these column indices belong.

\param col_indices (In/Out) Dense integer array filled with TRUE/FALSE values.  col_indices[j_box] is TRUE
       if the i_box equation depends on the j_box unknown.  On exit, all values should be set to FALSE.

\pre op must have been created by a call to create_linop.
\pre All column indices for all rows must have been inserted.

\post The graph associated with this DFT_LINOP object will be processed for optimal performance.
\post No further graph entries may be added after this function is called.
\post The attributes op->internal_index and op->external_index will be valid and usable.  

\return Returns 0 if no errors detected.

*/
int dft_fill_complete(DFT_LINOP * op);

/*! \fn int dft_insert_local_mat_row(DFT_LINOP * op, int i_mat, double * values)

\brief Insert values for i_mat row of matrix.

This function passes matrix values for the i_mat row to the DFT_LINOP object. The length of values
is assumed to be of compatible length with the number of columns on this processor as
determined by the call to fill_complete,
and in compatible ordering based on the arrays op->internal_index and op->external_index.

\param op (In) The address of the linear operator object that will be updated.

\param i_mat (In) Index of row in local "matrix" coordinate space to which these values belong.

\param values (In/Out) Dense integer array filled with matrix values.  values[j_mat] contains the
value of the i_mat row, j_mat column in local matrix coordinate space, compatible with
op->internal_index and op->external_index.  On exit, all values should be set to zero.

\pre op must have been created by a call to create_linop.
\pre All column indices for all rows must have been inserted and fill_complete called.

\post The values will be inserted in the matrix and nonzero locations in values array will be set to zero.  

\return Returns 0 if no errors detected.

*/
int dft_insert_local_mat_row(DFT_LINOP * op, int i_mat, double * values);


/*! \fn void dft_matvec(DFT_LINOP * op, double *x, double *y)

\brief Sparse matrix-vector multiply kernel.

\param A (In) A  pointer to a sparse matrix object.
\param x (In) Array of doubles of length compatible with the vec_space associated with the op.
\param y (Out) Array of doubles of length compatible with the vec_space associated with the op.

*/
void dft_matvec(DFT_LINOP * op, double *x, double *y);

#endif /* DFT_LINOP_H */
