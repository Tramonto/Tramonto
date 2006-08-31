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

#ifndef DFT_C2CPP_WRAPPERS_H
#define DFT_C2CPP_WRAPPERS_H
#include <mpi.h>

#include "dft_SolverManager.hpp"

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
