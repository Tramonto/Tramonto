
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

#ifdef EPETRA_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_C_wrappers.h"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#endif

#include "dft_c2cpp_wrappers.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*****************************************************/
  /**                  dft_EpetraComm                **/
  /***************************************************/

  DFT_OBJECT_PTR MANGLE(dft_epetrampicomm_create)(MPI_Comm * comm) {
    Epetra_Comm *comm_ = new Epetra_MpiComm(*comm);
    return((EPETRA_OBJECT_PTR ) comm_);
  }
  /*****************************************************/
  /**                  dft_SolverManager             **/
  /***************************************************/

  DFT_OBJECT_PTR MANGLE(dft_solvermanager_create)(DFT_INT numUnks, int* iunk_to_phys,
                        int* solverOptions, double* solverParams, DFT_OBJECT_REF comm) {
    dft_SolverManager * solvermanager_ = new dft_SolverManager(numUnks, iunk_to_phys, solverOptions,
		                                               solverParams, (MPI_Comm &) comm);
    return(solvermanager_);
  }

  void MANGLE(dft_solvermanager_destruct)(DFT_OBJECT_PTR solvermanagerptr) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanagerptr;
    solvermanager_->~dft_SolverManager();
  }

  int MANGLE(dft_solvermanager_setnodalrowmap)(DFT_OBJECT_PTR solvermanager, DFT_INT numgids, int * gids) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setNodalRowMap(numgids, gids));
  }

  int MANGLE(dft_solvermanager_setnodalcolmap)(DFT_OBJECT_PTR solvermanager, DFT_INT numgids, int * gids) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setNodalColMap(numgids, gids));
  }

  int MANGLE(dft_solvermanager_finalizeblockstructure)(DFT_OBJECT_PTR solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->finalizeBlockStructure());
  }

  int MANGLE(dft_solvermanager_initializeproblemvalues)(DFT_OBJECT_PTR solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->initializeProblemValues());
  }


  int MANGLE(dft_solvermanager_insertrhsvalue) (DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT inode, DFT_DOUBLE value) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertRhsValue(iunk, inode, value));
  }
  
  int MANGLE(dft_solvermanager_insertlhsvalue) (DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT inode, DFT_DOUBLE value) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertLhsValue(iunk, inode, value));
  }


  int MANGLE(dft_solvermanager_insertmatrixvalues) (DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT inode,
                                                    DFT_INT junk, DFT_INT numentries, double *values, int *indices) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertMatrixValues(iunk, inode, junk, numentries, values, indices));
  }

  int MANGLE(dft_solvermanager_insertonematrixvalue) (DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT inode,
                                                      DFT_INT junk, double value, int index) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertMatrixValues(iunk, inode, junk, 1, &value, &index));
  }
  
  int MANGLE(dft_solvermanager_finalizeproblemvalues)(DFT_OBJECT_PTR solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->finalizeProblemValues());
  }
  
  int MANGLE(dft_solvermanager_setblockmatrixreadonly)(DFT_OBJECT_PTR solvermanager, DFT_INT iunk, DFT_INT junk, DFT_INT readOnly) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    bool readonly_ = !(DFT_DEREF(readOnly==0));
    return(solvermanager_->setBlockMatrixReadOnly(iunk, junk, readonly_));
  }

  int MANGLE(dft_solvermanager_setlhs)(DFT_OBJECT_PTR solvermanager, double** x) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setLhs((const double**) x));
  }

  int MANGLE(dft_solvermanager_setrhs)(DFT_OBJECT_PTR solvermanager, double** x) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setRhs((const double**) x));
  }

  int MANGLE(dft_solvermanager_getlhs)(DFT_OBJECT_PTR solvermanager, double** x) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->getLhs(x));
  }

  int MANGLE(dft_solvermanager_getrhs)(DFT_OBJECT_PTR solvermanager, double** x) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->getRhs(x));
  }

  int MANGLE(dft_solvermanager_setupsolver)(DFT_OBJECT_PTR solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setupSolver());
  }

  int MANGLE(dft_solvermanager_solve)(DFT_OBJECT_PTR solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->solve());
  }

  int MANGLE(dft_solvermanager_importr2c)(DFT_OBJECT_PTR solvermanager, double** x, double **b) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->importR2C((const double**) x,b));
  }

  int MANGLE(dft_solvermanager_importnodalr2c)(DFT_OBJECT_PTR solvermanager, double* x, double *b) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->importR2C((const double*) x,b));
  }

#ifdef __cplusplus
}
#endif
