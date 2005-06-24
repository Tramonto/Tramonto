
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER


#include "dft_c2cpp_wrappers.h"
#include "dft_BasicLinProbMgr.hpp"

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_BasicLinProbMgr             **/
  /***************************************************/

  void * dft_linprobmgr_create(int numUnks,
                        int* solverOptions, double* solverParams, MPI_Comm comm) {
    dft_BasicLinProbMgr * linprobmgr_ = new dft_BasicLinProbMgr(numUnks, solverOptions,
		                                               solverParams, comm);
    return((void *)linprobmgr_);
  }

  void dft_linprobmgr_destruct(void * linprobmgr) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    delete linprobmgr_;
  }

  int dft_linprobmgr_setnodalrowmap(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->setNodalRowMap(numgids, gids));
  }

  int dft_linprobmgr_setnodalcolmap(void * linprobmgr, int numgids, int * gids) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->setNodalColMap(numgids, gids));
  }

  int dft_linprobmgr_finalizeblockstructure(void * linprobmgr) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->finalizeBlockStructure());
  }

  int dft_linprobmgr_initializeproblemvalues(void * linprobmgr) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->initializeProblemValues());
  }


  int dft_linprobmgr_insertrhsvalue (void * linprobmgr, int iunk, int inode, double value) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->insertRhsValue(iunk, inode, value));
  }

  int dft_linprobmgr_insertonematrixvalue (void * linprobmgr, int iunk, int ownednode,
                                                      int junk, int boxnode, double value) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->insertMatrixValue(iunk, ownednode, junk, boxnode, value));
  }
  
  int dft_linprobmgr_insertmultinodematrixvalues (void * linprobmgr, int iunk, int ownednode,
                                                    int junk, int * boxnodeindices, double *values, int numEntries) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->insertMatrixValues(iunk, ownednode, junk, boxnodeindices, values, numEntries));
  }
  
  int dft_linprobmgr_insertmultiphysicsmatrixvalues (void * linprobmgr, int iunk, int ownednode,
                                                    int *junkindices, int boxnode, double *values, int numEntries) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->insertMatrixValues(iunk, ownednode, junkindices, boxnode, values, numEntries));
  }

  int dft_linprobmgr_finalizeproblemvalues(void * linprobmgr) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->finalizeProblemValues());
  }
  
  int dft_linprobmgr_setblockmatrixreadonly(void * linprobmgr, int iunk, int junk, int readOnly) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    bool readonly_ = !(readOnly==0);
    return(linprobmgr_->setBlockMatrixReadOnly(iunk, junk, readonly_));
  }

  int dft_linprobmgr_setrhs(void * linprobmgr, double** x) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->setRhs((const double**) x));
  }

  int dft_linprobmgr_getlhs(void * linprobmgr, double** x) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->getLhs(x));
  }

  int dft_linprobmgr_getrhs(void * linprobmgr, double** x) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->getRhs(x));
  }

  int dft_linprobmgr_writeMatrix(void * linprobmgr, char * filename, char * matrixName, char * matrixDescription) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->writeMatrix(filename, matrixName, matrixDescription));
  }

  int dft_linprobmgr_setupsolver(void * linprobmgr) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->setupSolver());
  }

  int dft_linprobmgr_solve(void * linprobmgr) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->solve());
  }

  int dft_linprobmgr_applymatrix(void * linprobmgr, double**x, double** b) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->applyMatrix((const double**) x,b));
  }

  int dft_linprobmgr_importr2c(void * linprobmgr, double** x, double **b) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->importR2C((const double**) x,b));
  }

  int dft_linprobmgr_importnodalr2c(void * linprobmgr, double* x, double *b) {
    dft_BasicLinProbMgr * linprobmgr_ = (dft_BasicLinProbMgr *) linprobmgr;
    return(linprobmgr_->importR2C((const double*) x,b));
  }

#ifdef __cplusplus
}
#endif
