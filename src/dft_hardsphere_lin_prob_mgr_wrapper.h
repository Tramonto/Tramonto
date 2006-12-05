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

#ifndef DFT_HARDSPHERE_LIN_PROB_MGR_WRAPPER_H
#define DFT_HARDSPHERE_LIN_PROB_MGR_WRAPPER_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_HardSphereLinProbMgr             **/
  /***************************************************/

  void * dft_hardsphere_lin_prob_mgr_create(int numUnks,
				      int* solverOptions, double* solverParams, MPI_Comm comm);

  void * dft_hardsphere_lin_prob_mgr_create_debug(int numUnks,
					    int* solverOptions, double* solverParams, MPI_Comm comm);

  void dft_hardsphere_lin_prob_mgr_destruct(void * linprobmgr);

  int dft_hardsphere_lin_prob_mgr_setindnonlocalequationids(void * linprobmgr, int numgids, int * gids);
  
  int dft_hardsphere_lin_prob_mgr_setdepnonlocalequationids(void * linprobmgr, int numgids, int * gids);
  
  int dft_hardsphere_lin_prob_mgr_setdensityequationids(void * linprobmgr, int numgids, int * gids);

  int dft_hardsphere_lin_prob_mgr_seta22blockisdiagona(void * linprobmgr, int isa22diagonal);

#ifdef __cplusplus
}
#endif

#endif /* DFT_HARDSPHERE_LIN_PROB_MGR_WRAPPER_H */
