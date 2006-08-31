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

#ifndef DFT_LINOP_H
#define DFT_LINOP_H

#include "dft_vecspace.h"

/*! \file dft_linop.h
\brief File defining struct that contains all information about a linear operator object.

The file must be included by each Tramonto source file that needs access
to linear operators.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \struct DFT_LINOP_STRUCT dft_linop.h

\brief Contains all information about a linear operator object.

*/
struct DFT_LINOP_STRUCT {

  DFT_VECSPACE * vecspace;    /*!< Selector for linear algebra library */
  int first_time_graph_insert;   /*!< Initialized to TRUE, set to FALSE after the first call to insert_global_graph_row */
  int * internal_index; /*!< After fill_complete is called, contains the operator ordering for GIDs owned by this processor */
  int * external_index; /*!< After fill_complete is called, contains the operator ordering for GIDs not owned by this processor but present as column indices on this processor */
  int nonzeros;     /*!< Number of nonzeros in linear operator */
  void * aux_info;  /*!< Generic pointer that can be used by the implementing library for addional data */
  
};

/*! \typedef DFT_LINOP

\brief Convenient expression for struct DFT_LINOP_STRUCT.

*/
typedef struct DFT_LINOP_STRUCT DFT_LINOP;

/*! \fn void dft_create_linop(DFT_VECSPACE * vecspace, DFT_LINOP ** linop)

\brief Create a vector space of size N_loc on each processor with global IDs from update.

\param vecspace (In) Vector space containing information to create operator.
\param linop (Out) The address of the linear operator.

\pre vecspace must be properly defined.
\pre linop must be an uninitialized pointer.  If it is already initialized, there will be 
a memory leak.

\post The vecspace pointer will be copied, but the associated vector space data WILL NOT be copied.  It is the
responsibility of the calling program to maintain the integrity of the associated vector space object.

*/

void dft_create_linop(DFT_VECSPACE * vecspace, DFT_LINOP ** linop);

/*! \fn void dft_destroy_linop(DFT_LINOP ** linop)

\brief Destroy a linear operator and delete all associated memory.

\param linop (In/Out) The address of the linear operator object that will be destroyed.

\pre linop must have been created by a call to dft_create_linop.

\post All memory associated with this object will be deleted and the pointer will be set to NULL.  
\post The associated vecspace object will NOT be deleted.  This is the responsibility of the calling
program.

*/

void dft_destroy_linop(DFT_LINOP ** linop);

/*! \fn int dft_insert_global_graph_row(DFT_LINOP * linop, int i_mat, int * indices, int nonzeros_in_row)

\brief Insert box column indices for i_box row into i_mat row of matrix graph.

This function passes column indices for the graph of the matrix to the DFT_LINOP object.  The column 
indices are in box coordinate space.

\param linop (In) The address of the linear operator object that will be updated.

\param i_mat (In) Index of row in local "matrix" coordinate space to which these column indices belong.

\param indices (In) Array of column indices that are have nonzero entries in this row.

\pre linop must have been created by a call to create_linop.

\post The graph associated with this DFT_LINOP object will be updated with the passed-in column entries.
\post This row of the matrix may not receive additional column indices and the number of column indices is fixed.

\return Returns 0 if no errors detected.

*/
int dft_insert_global_graph_row(DFT_LINOP * linop, int i_mat, int * indices, int nonzeros_in_row);

/*! \fn int dft_fill_complete(DFT_LINOP * linop)

\brief Insert box column indices for i_box row into i_mat row of matrix graph.

This function passes column indices for the graph of the matrix to the DFT_LINOP object.  The column 
indices are in box coordinate space.

\param linop (In) The address of the linear operator object that will be updated.

\param i_box (In) Index of row in "box" coordinate space to which these column indices belong.

\param i_mat (In) Index of row in local "matrix" coordinate space to which these column indices belong.

\param col_indices (In/Out) Dense integer array filled with TRUE/FALSE values.  col_indices[j_box] is TRUE
       if the i_box equation depends on the j_box unknown.  On exit, all values should be set to FALSE.

\pre linop must have been created by a call to create_linop.
\pre All column indices for all rows must have been inserted.

\post The graph associated with this DFT_LINOP object will be processed for optimal performance.
\post No further graph entries may be added after this function is called.
\post The attributes linop->internal_index and linop->external_index will be valid and usable.  

\return Returns 0 if no errors detected.

*/
int dft_fill_complete(DFT_LINOP * linop);

/*! \fn int dft_insert_local_mat_row(DFT_LINOP * linop, int i_mat, double * values)

\brief Insert values for i_mat row of matrix.

This function passes matrix values for the i_mat row to the DFT_LINOP object. The length of values
is assumed to be of compatible length with the number of columns on this processor as
determined by the call to fill_complete,
and in compatible ordering based on the arrays linop->internal_index and linop->external_index.

\param linop (In) The address of the linear operator object that will be updated.

\param i_mat (In) Index of row in local "matrix" coordinate space to which these values belong.

\param values (In/Out) Dense integer array filled with matrix values.  values[j_mat] contains the
value of the i_mat row, j_mat column in local matrix coordinate space, compatible with
op->internal_index and linop->external_index.  On exit, all values should be set to zero.

\pre linop must have been created by a call to create_linop.
\pre All column indices for all rows must have been inserted and fill_complete called.

\post The values will be inserted in the matrix and nonzero locations in values array will be set to zero.  

\return Returns 0 if no errors detected.

*/
int dft_insert_local_mat_row(DFT_LINOP * linop, int i_mat, double * values);


/*! \fn void dft_matvec(DFT_LINOP * linop, double *x, double *y)

\brief Sparse matrix-vector multiply kernel.

\param linop (In) A  pointer to a sparse matrix object.
\param x (In) Array of doubles of length compatible with the vecspace associated with the linop.
\param y (Out) Array of doubles of length compatible with the vecspace associated with the linop.

*/
void dft_matvec(DFT_LINOP * linop, double *x, double *y);

/*! \fn void dft_create_range_vec(DFT_LINOP * linop, double ** y)

\brief Create a vector compatible with the range of the linear operator, eg, a vector y such that y = A*x.

\param linop (In) A  pointer to a sparse matrix object.
\param y (Out) Array of doubles of length compatible with the range associated with the linop.  For Aztec and
Epetra, this vector is typically the same size on each processor as the number of rows kept on the processor.
*/
void dft_create_range_vec(DFT_LINOP * linop, double **y);

/*! \fn void dft_create_domain_vec(DFT_LINOP * linop, double ** x)

\brief Create a vector compatible with the domain of the linear operator, eg, a vector x such that y = A*x.

\param linop (In) A  pointer to a sparse matrix object.
\param x (Out) Array of doubles of length compatible with the domain associated with the linop.  For Aztec,
this vector is typically longer than the number of rows on the processor.  The extra length is used to store
ghost node values for sparse matrix vector multiplication.  For Epetra this vector is typically the same length
as the number of rows on the calling processor.
*/
void dft_create_domain_vec(DFT_LINOP * linop, double **x);

#endif /* DFT_LINOP_H */
