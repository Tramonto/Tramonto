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


#include "dft_poly_lin_prob_mgr_wrapper.h"
#include "dft_PolyLinProbMgr.hpp"

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_BasicLinProbMgr             **/
  /***************************************************/

  void * dft_poly_lin_prob_mgr_create(int numUnks,
                        int* solverOptions, double* solverParams, MPI_Comm comm) {
    dft_PolyLinProbMgr * tmp = new dft_PolyLinProbMgr(numUnks, solverOptions,
		                                               solverParams, comm);
    dft_BasicLinProbMgr * linprobmgr_ = dynamic_cast<dft_BasicLinProbMgr *>(tmp);
    return((void *)linprobmgr_);
  }

  void * dft_poly_lin_prob_mgr_create_debug(int numUnks,
                        int* solverOptions, double* solverParams, MPI_Comm comm) {
    dft_PolyLinProbMgr * tmp = new dft_PolyLinProbMgr(numUnks, solverOptions,
						      solverParams, comm, true);
    dft_BasicLinProbMgr * linprobmgr_ = dynamic_cast<dft_BasicLinProbMgr *>(tmp);
    return((void *)linprobmgr_);
  }

  void dft_poly_lin_prob_mgr_destruct(void * linprobmgr) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_PolyLinProbMgr * linprobmgr_ = dynamic_cast<dft_PolyLinProbMgr *>(tmp);
    delete linprobmgr_;
  }

  int dft_poly_lin_prob_mgr_setgequationids(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_PolyLinProbMgr * linprobmgr_  = dynamic_cast<dft_PolyLinProbMgr *>(tmp);
    return(linprobmgr_->setGEquationIDs(numgids, gids));
  }

  int dft_poly_lin_prob_mgr_setginvequationids(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_PolyLinProbMgr * linprobmgr_  = dynamic_cast<dft_PolyLinProbMgr *>(tmp);
    return(linprobmgr_->setGInvEquationIDs(numgids, gids));
  }

  int dft_poly_lin_prob_mgr_setcmsequationids(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_PolyLinProbMgr * linprobmgr_  = dynamic_cast<dft_PolyLinProbMgr *>(tmp);
    return(linprobmgr_->setCmsEquationIDs(numgids, gids));
  }

  int dft_poly_lin_prob_mgr_setdensityequationids(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_PolyLinProbMgr * linprobmgr_  = dynamic_cast<dft_PolyLinProbMgr *>(tmp);
    return(linprobmgr_->setDensityEquationIDs(numgids, gids));
  }

  int dft_poly_lin_prob_mgr_setfieldondensityislinear(void * linprobmgr, int isLinear) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_PolyLinProbMgr * linprobmgr_  = dynamic_cast<dft_PolyLinProbMgr *>(tmp);
    bool isLinear1 = (isLinear!=0);
    return(linprobmgr_->setFieldOnDensityIsLinear(isLinear1));
  }

#ifdef __cplusplus
}
#endif
