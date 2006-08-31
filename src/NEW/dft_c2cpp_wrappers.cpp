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


#include "dft_c2cpp_wrappers.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_SolverManager             **/
  /***************************************************/

  void * dft_solvermanager_create(int numUnks, int* Unk2Phys,
                        int* solverOptions, double* solverParams, MPI_Comm comm) {
    dft_SolverManager * solvermanager_ = new dft_SolverManager(numUnks, Unk2Phys, solverOptions,
		                                               solverParams, comm);
    return((void *)solvermanager_);
  }

  void dft_solvermanager_destruct(void * solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    delete solvermanager_;
  }

  int dft_solvermanager_setnodalrowmap(void * solvermanager, int numgids, int * gids) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setNodalRowMap(numgids, gids));
  }

  int dft_solvermanager_setnodalcolmap(void * solvermanager, int numgids, int * gids) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setNodalColMap(numgids, gids));
  }

  int dft_solvermanager_finalizeblockstructure(void * solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->finalizeBlockStructure());
  }

  int dft_solvermanager_initializeproblemvalues(void * solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->initializeProblemValues());
  }


  int dft_solvermanager_insertrhsvalue (void * solvermanager, int iunk, int inode, double value) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertRhsValue(iunk, inode, value));
  }

  int dft_solvermanager_insertonematrixvalue (void * solvermanager, int iunk, int ownednode,
                                                      int junk, int boxnode, double value) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertMatrixValue(iunk, ownednode, junk, boxnode, value));
  }
  
  int dft_solvermanager_insertmultinodematrixvalues (void * solvermanager, int iunk, int ownednode,
                                                    int junk, int * boxnodeindices, double *values, int numEntries) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertMatrixValues(iunk, ownednode, junk, boxnodeindices, values, numEntries));
  }
  
  int dft_solvermanager_insertmultiphysicsmatrixvalues (void * solvermanager, int iunk, int ownednode,
                                                    int *junkindices, int boxnode, double *values, int numEntries) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertMatrixValues(iunk, ownednode, junkindices, boxnode, values, numEntries));
  }

  int dft_solvermanager_finalizeproblemvalues(void * solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->finalizeProblemValues());
  }
  
  int dft_solvermanager_setblockmatrixreadonly(void * solvermanager, int iunk, int junk, int readOnly) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    bool readonly_ = !(readOnly==0);
    return(solvermanager_->setBlockMatrixReadOnly(iunk, junk, readonly_));
  }

  int dft_solvermanager_setrhs(void * solvermanager, double** x) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setRhs((const double**) x));
  }

  int dft_solvermanager_getlhs(void * solvermanager, double** x) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->getLhs(x));
  }

  int dft_solvermanager_getrhs(void * solvermanager, double** x) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->getRhs(x));
  }

  int dft_solvermanager_setupsolver(void * solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setupSolver());
  }

  int dft_solvermanager_solve(void * solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->solve());
  }

  int dft_solvermanager_importr2c(void * solvermanager, double** x, double **b) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->importR2C((const double**) x,b));
  }

  int dft_solvermanager_importnodalr2c(void * solvermanager, double* x, double *b) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->importR2C((const double*) x,b));
  }

#ifdef __cplusplus
}
#endif
