
//@HEADER
/*
************************************************************************

Epetra: Linear Algebra Services Package 
Copyright (2001) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Michael A. Heroux (maherou@sandia.gov) 

************************************************************************
*/
//@HEADER

#ifndef DFT_C2CPP_WRAPPERS_H
#define DFT_C2CPP_WRAPPERS_H

#ifdef DFT_FORTRAN

typedef double * DFT_DOUBLE;
typedef int    * DFT_INT;
#define DFT_DEREF(a) *a

#ifdef DFT_ADDRESS64BIT

typedef long int DFT_OBJECT_PTR;
typedef long int & DFT_OBJECT_REF;

#else

typedef int DFT_OBJECT_PTR;
typedef int & DFT_OBJECT_REF;

#endif
#else

/* These typedefs act as new types for the Epetra C interface */
;
typedef double DFT_DOUBLE;
typedef int    DFT_INT;

#define DFT_DEREF(a) a

typedef void * DFT_OBJECT_PTR;
typedef void * DFT_OBJECT_REF;


#endif
 
#ifdef DFT_FORTRAN
#if defined(DFT_HAVE_NO_FORTRAN_UNDERSCORE)
#define MANGLE(x) x
#else
#define MANGLE(x) x ## __
#endif
#else
#define MANGLE(x) x
#endif


#include "dft_SolverManager.hpp"

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_EpetraComm                **/
  /***************************************************/

  DFT_OBJECT_PTR MANGLE(dft_epetrampicomm_create)(MPI_Comm * comm);
  /*****************************************************/
  /**                  dft_SolverManager             **/
  /***************************************************/

  DFT_OBJECT_PTR MANGLE(dft_solvermanager_create)(DFT_INT numUnks, int* iunk_to_phys,
		        int* solverOptions, double* solverParams, DFT_OBJECT_REF comm);

  void MANGLE(dft_solvermanager_destruct)(DFT_OBJECT_PTR solvermanagerptr);

  int MANGLE(dft_solvermanager_setnodalrowmap)(DFT_OBJECT_PTR solvermanager, DFT_INT numgids, int * gids);

  int MANGLE(dft_solvermanager_setnodalcolmap)(DFT_OBJECT_PTR solvermanager, DFT_INT numgids, int * gids);

  int MANGLE(dft_solvermanager_finalizeblockstructure)(DFT_OBJECT_PTR solvermanager);

  int MANGLE(dft_solvermanager_initializeproblemvalues)(DFT_OBJECT_PTR solvermanager);
  
  int MANGLE(dft_solvermanager_insertrhsvalue)(DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT inode, DFT_DOUBLE value);

  int MANGLE(dft_solvermanager_insertlhsvalue)(DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT ibox, DFT_DOUBLE value);
  
  int MANGLE(dft_solvermanager_insertmatrixvalues)(DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT inode,
		                                   DFT_INT junk, DFT_INT numentries, double *values, int *indices);

  int MANGLE(dft_solvermanager_insertonematrixvalue)(DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT inode,
                                                     DFT_INT junk, double value, int index);
  
  int MANGLE(dft_solvermanager_finalizeproblemvalues)(DFT_OBJECT_PTR solvermanager);
  
  int MANGLE(dft_solvermanager_setblockmatrixreadonly)(DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT junk, DFT_INT readOnly);

  int MANGLE(dft_solvermanager_setlhs)(DFT_OBJECT_PTR solvermanager, double** x);
  
  int MANGLE(dft_solvermanager_setrhs)(DFT_OBJECT_PTR solvermanager, double** b);

  int MANGLE(dft_solvermanager_getlhs)(DFT_OBJECT_PTR solvermanager, double** x);

  int MANGLE(dft_solvermanager_getrhs)(DFT_OBJECT_PTR solvermanager, double** b);

  int MANGLE(dft_solvermanager_setupsolver)(DFT_OBJECT_PTR solvermanager);

  int MANGLE(dft_solvermanager_solve)(DFT_OBJECT_PTR solvermanager);

  int MANGLE(dft_solvermanager_applymatrix)(DFT_OBJECT_PTR solvermanager, double**x, double** b);

  int MANGLE(dft_solvermanager_importr2c)(DFT_OBJECT_PTR solvermanager, double**x, double** b);

  int MANGLE(dft_solvermanager_importnodalr2c)(DFT_OBJECT_PTR solvermanager, double*x, double* b);


#ifdef __cplusplus
}
#endif

#endif /* DFT_C2CPP_WRAPPERS_H */
