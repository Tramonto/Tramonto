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

#ifndef DFT_POLY_LIN_PROB_MGR_WRAPPER_H
#define DFT_POLY_LIN_PROB_MGR_WRAPPER_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_PolyLinProbMgr             **/
  /***************************************************/

  void * dft_poly_lin_prob_mgr_create(int numUnks,
				      int* solverOptions, double* solverParams, MPI_Comm comm);

  void * dft_poly_lin_prob_mgr_create_debug(int numUnks,
					    int* solverOptions, double* solverParams, MPI_Comm comm);

  void dft_poly_lin_prob_mgr_destruct(void * linprobmgr);

  int dft_poly_lin_prob_mgr_setgequationids(void * linprobmgr, int numgids, int * gids);
  
  int dft_poly_lin_prob_mgr_setginvequationids(void * linprobmgr, int numgids, int * gids);
  
  int dft_poly_lin_prob_mgr_setcmsequationids(void * linprobmgr, int numgids, int * gids);
  
  int dft_poly_lin_prob_mgr_setdensityequationids(void * linprobmgr, int numgids, int * gids);
  int dft_poly_lin_prob_mgr_setfieldondensityislinear(void * linprobmgr, int isLinear);

#ifdef __cplusplus
}
#endif

#endif /* DFT_POLY_LIN_PROB_MGR_WRAPPER_H */
