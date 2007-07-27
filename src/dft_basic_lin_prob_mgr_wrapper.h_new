
/*
//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER
*/

#ifndef DFT_BASIC_LIN_PROB_MGR_WRAPPER_H
#define DFT_BASIC_LIN_PROB_MGR_WRAPPER_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_Linprobmgr             **/
  /***************************************************/

  /*  void * dft_basic_lin_prob_mgr_create(int numUnks,
      int* solverOptions, double* solverParams, MPI_Comm comm);*/

  void * dft_basic_lin_prob_mgr_create(int numUnks, void * Parameterlist_list, MPI_Comm comm);

  void dft_linprobmgr_destruct(void * linprobmgrptr);

  int dft_linprobmgr_setnodalrowmap(void * linprobmgr, int numgids, int * gids);

  int dft_linprobmgr_setnodalcolmap(void * linprobmgr, int numgids, int * gids);

  int dft_linprobmgr_setcoarsenednodeslist(void * linprobmgr, int numgids, int * gids);

  int dft_linprobmgr_finalizeblockstructure(void * linprobmgr);

  int dft_linprobmgr_initializeproblemvalues(void * linprobmgr);
  
  int dft_linprobmgr_insertrhsvalue(void * linprobmgr, int iunk, int inode, double value);

  int dft_linprobmgr_insertonematrixvalue(void * linprobmgr, int iunk, int ownednode,
                                                     int junk, int boxnode, double value);

  int dft_linprobmgr_insertmultinodematrixvalues(void * linprobmgr, int iunk, int ownednode,
							    int junk, int *boxnodeindices, double *values, int numentries);

  int dft_linprobmgr_insertmultiphysicsmatrixvalues(void * linprobmgr, int iunk, int ownednode,
		                                   int * junkindices, int boxnode, double *values, int numentries);
  
  int dft_linprobmgr_finalizeproblemvalues(void * linprobmgr);
  
  double dft_linprobmgr_getmatrixvalue(void * linprobmgr, int iunk, int ownednode,
                                                     int junk, int boxnode);

  int dft_linprobmgr_setrhs(void * linprobmgr, double** b);

  int dft_linprobmgr_getlhs(void * linprobmgr, double** x);

  int dft_linprobmgr_getrhs(void * linprobmgr, double** b);

  int dft_linprobmgr_writeMatrix(void * linprobmgr, char * filename, char * matrixName, char * matrixDescription);

  int dft_linprobmgr_setupsolver(void * linprobmgr);

  int dft_linprobmgr_solve(void * linprobmgr);

  int dft_linprobmgr_applymatrix(void * linprobmgr, double**x, double** b);

  int dft_linprobmgr_importr2c(void * linprobmgr, double**x, double** b);

  int dft_linprobmgr_importnodalr2c(void * linprobmgr, double*x, double* b);


#ifdef __cplusplus
}
#endif

#endif /* DFT_BASIC_LIN_PROB_MGR_WRAPPER_H */
