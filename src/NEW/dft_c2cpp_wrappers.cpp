
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

  DFT_OBJECT_PTR MANGLE(dft_solvermanager_create)(DFT_INT numblocks, DFT_OBJET_REF comm) {
    Epetra_Comm& comm_ = *(Epetra_Comm *) comm;
    Epetra_Map *map = new Epetra_Map(EPETRA_DEREF(numGlobalElements), EPETRA_DEREF(indexBase), comm_);
    return((EPETRA_OBJECT_PTR ) map);
  }

  int MANGLE(dft_solvermanager_setrowmap)(DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT numgids, int * gids) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setRowMap(i, numgids, gids));
  }

  int MANGLE(dft_solvermanager_setcolmap)(DFT_OBJECT_REF solvermanager, DFT_INT j, DFT_INT numgids, int * gids) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->setColMap(j, numgids, gids));
  }

  int MANGLE(dft_solvermanager_finalizeblockstructure)(DFT_OBJECT_REF solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->finalizeBlockStructure());
  }

  int  MANGLE(dft_solvermanager_insertgraphindices) (DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT j, DFT_INT localrow, DFT_INT numentries, int *indices) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertGraphIndices(i, j, localrow, numentries, indices));
  }

  int MANGLE(dft_solvermanager_finalizegraphstructure)(DFT_OBJECT_REF solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->finalizeGraphStructure());
  }

  int MANGLE(dft_solvermanager_initializeproblemvalues)(DFT_OBJECT_REF solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->initializeProblemValues());
  }

  
  int MANGLE(dft_solvermanager_insertmatrixvalues) (DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT j, DFT_INT localrow, DFT_INT numentries, double *values, int *indices) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertMatrixValues(i, j, localrow, numentries, values, indices));
  }

  
  int MANGLE(dft_solvermanager_insertlhsvalues) (DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT localentry, DFT_DOUBLE value) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertLhsValues(i, localentry, value));
  }

  
  int MANGLE(dft_solvermanager_insertrhsvalues) (DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT localentry, DFT_DOUBLE value) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->insertRhsValues(i, localentry, value));
  }
  
  int MANGLE(dft_solvermanager_finalizeproblemvalues)(DFT_OBJECT_REF solvermanager) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    return(solvermanager_->finalizeProblemValues());
  }

  
  int MANGLE(dft_solvermanager_setblockmatrixreadonly)(DFT_OBJECT_REF solvermanager, DFT_INT i, DFT_INT j, DFT_INT readOnly) {
    dft_SolverManager * solvermanager_ = (dft_SolverManager *) solvermanager;
    bool readonly_ = !(DFT_DEREF(readOnly==0));
    return(solvermanager_->setblockmatrixreadonly(i, j, readonly_));
  }

#ifdef __cplusplus
}
#endif
