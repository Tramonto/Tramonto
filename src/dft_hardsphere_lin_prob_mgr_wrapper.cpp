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


#include "dft_hardsphere_lin_prob_mgr_wrapper.h"
#include "dft_HardSphereLinProbMgr.hpp"

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_HardSphereLinProbMgr             **/
  /***************************************************/
  void * dft_hardsphere_lin_prob_mgr_create(int numUnks,
                        int* solverOptions, double* solverParams, MPI_Comm comm) {
    dft_HardSphereLinProbMgr * tmp = new dft_HardSphereLinProbMgr(numUnks, solverOptions,
		                                               solverParams, comm);
    dft_BasicLinProbMgr * linprobmgr_ = dynamic_cast<dft_BasicLinProbMgr *>(tmp);
    return((void *)linprobmgr_);
  }

  void * dft_hardsphere_lin_prob_mgr_create_debug(int numUnks,
                        int* solverOptions, double* solverParams, MPI_Comm comm) {
    dft_HardSphereLinProbMgr * tmp = new dft_HardSphereLinProbMgr(numUnks, solverOptions,
						      solverParams, comm, true);
    dft_BasicLinProbMgr * linprobmgr_ = dynamic_cast<dft_BasicLinProbMgr *>(tmp);
    return((void *)linprobmgr_);
  }

  void dft_hardsphere_lin_prob_mgr_destruct(void * linprobmgr) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_HardSphereLinProbMgr * linprobmgr_ = dynamic_cast<dft_HardSphereLinProbMgr *>(tmp);
    delete linprobmgr_;
  }

  int dft_hardsphere_lin_prob_mgr_setindnonlocalequationids(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_HardSphereLinProbMgr * linprobmgr_  = dynamic_cast<dft_HardSphereLinProbMgr *>(tmp);
    return(linprobmgr_->setIndNonLocalEquationIDs(numgids, gids));
  }

  int dft_hardsphere_lin_prob_mgr_setdepnonlocalequationids(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_HardSphereLinProbMgr * linprobmgr_  = dynamic_cast<dft_HardSphereLinProbMgr *>(tmp);
    return(linprobmgr_->setDepNonLocalEquationIDs(numgids, gids));
  }

  int dft_hardsphere_lin_prob_mgr_setdensityequationids(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_HardSphereLinProbMgr * linprobmgr_  = dynamic_cast<dft_HardSphereLinProbMgr *>(tmp);
    return(linprobmgr_->setDensityEquationIDs(numgids, gids));
  }

  int dft_hardsphere_lin_prob_mgr_seta22blockisdiagonal(void * linprobmgr, int isa22diagonal) {
    dft_BasicLinProbMgr * tmp = (dft_BasicLinProbMgr *) linprobmgr;
    dft_HardSphereLinProbMgr * linprobmgr_  = dynamic_cast<dft_HardSphereLinProbMgr *>(tmp);
    bool isA22Diagonal = (isa22diagonal!=0);
    return(linprobmgr_->setA22BlockIsDiagonal(isA22Diagonal));
  }

#ifdef __cplusplus
}
#endif
