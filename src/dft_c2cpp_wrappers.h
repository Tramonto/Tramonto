
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
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_SolverManager             **/
  /***************************************************/

  void * dft_solvermanager_create(int numUnks, int* Unk2Phys,
		        int* solverOptions, double* solverParams, MPI_Comm comm);

  void dft_solvermanager_destruct(void * solvermanagerptr);

  int dft_solvermanager_setnodalrowmap(void * solvermanager, int numgids, int * gids);

  int dft_solvermanager_setnodalcolmap(void * solvermanager, int numgids, int * gids);

  int dft_solvermanager_finalizeblockstructure(void * solvermanager);

  int dft_solvermanager_initializeproblemvalues(void * solvermanager);
  
  int dft_solvermanager_insertrhsvalue(void * solvermanager, int iunk, int inode, double value);

  int dft_solvermanager_insertonematrixvalue(void * solvermanager, int iunk, int ownednode,
                                                     int junk, int boxnode, double value);

  int dft_solvermanager_insertmultinodematrixvalues(void * solvermanager, int iunk, int ownednode,
							    int junk, int *boxnodeindices, double *values, int numentries);

  int dft_solvermanager_insertmultiphysicsmatrixvalues(void * solvermanager, int iunk, int ownednode,
		                                   int * junkindices, int boxnode, double *values, int numentries);
  
  int dft_solvermanager_finalizeproblemvalues(void * solvermanager);
  
  int dft_solvermanager_setblockmatrixreadonly(void * solvermanager, int iunk, int junk, int readOnly);
  
  int dft_solvermanager_setrhs(void * solvermanager, double** b);

  int dft_solvermanager_getlhs(void * solvermanager, double** x);

  int dft_solvermanager_getrhs(void * solvermanager, double** b);

  int dft_solvermanager_setupsolver(void * solvermanager);

  int dft_solvermanager_solve(void * solvermanager);

  int dft_solvermanager_applymatrix(void * solvermanager, double**x, double** b);

  int dft_solvermanager_importr2c(void * solvermanager, double**x, double** b);

  int dft_solvermanager_importnodalr2c(void * solvermanager, double*x, double* b);


#ifdef __cplusplus
}
#endif

#endif /* DFT_C2CPP_WRAPPERS_H */
