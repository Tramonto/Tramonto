
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

  DFT_OBJECT_PTR MANGLE(dft_solvermanager_create)(DFT_INT numblocks, DFT_OBJET_REF comm);

  int MANGLE(dft_solvermanager_setrowmap)(DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT numgids, int * gids);

  int MANGLE(dft_solvermanager_setcolmap)(DFT_OBJECT_REF solvermanager, DFT_INT j, DFT_INT numgids, int * gids);

  int MANGLE(dft_solvermanager_finalizeblockstructure)(DFT_OBJECT_REF solvermanager);

  int  MANGLE(dft_solvermanager_insertgraphindices) (DFT_INT i, DFT_INT j, DFT_INT localrow, DFT_INT numentries, int *indices);

  int MANGLE(dft_solvermanager_finalizegraphstructure)(DFT_OBJECT_REF solvermanager);
  
  int MANGLE(dft_solvermanager_initializeproblemvalues)(DFT_OBJECT_REF solvermanager);
  
  int MANGLE(dft_solvermanager_insertmatrixvalues) (DFT_INT i, DFT_INT j, DFT_INT localrow, DFT_INT numentries, double *values, int *indices);
  
  int MANGLE(dft_solvermanager_insertlhsvalues) (DFT_INT i, DFT_INT localentry, DFT_DOUBLE value);
  
  int MANGLE(dft_solvermanager_insertrhsvalues) (DFT_INT i, DFT_INT localentry, DFT_DOUBLE value);
  
  int MANGLE(dft_solvermanager_finalizeproblemvalues)(DFT_OBJECT_REF solvermanager);
  
  int MANGLE(dft_solvermanager_setblockmatrixreadonly)(DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT j, DFT_INT readOnly);

#ifdef __cplusplus
}
#endif

#endif /* DFT_C2CPP_WRAPPERS_H */
