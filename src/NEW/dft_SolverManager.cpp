/*@HEADER
// ***********************************************************************
// 
//                Tramonto: Molecular Theories Modeling Code
//                 Copyright (2004) Sandia Corporation
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
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// Questions? Contact Laura J.D. Frink (ljfrink@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/

#include "dft_SolverManager.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"


//=============================================================================
dft_SolverManager::dft_SolverManager(int numBlocks, const Epetra_Comm & comm) 
  : numBlocks_(numBlocks),
    comm_(comm),
    blockMatrix_(0),
    blockGraph_(0),
    blockLhs_(0),
    blockRhs_(0),
    rowMaps_(0),
    colMaps_(0),
    isBlockStructureSet_(false),
    isGraphStructureSet_(false) {
  // Allocate blockMatrix, blockGraph, blockLhs, blockRhs and map arrays
  blockMatrix_ = new Epetra_CrsMatrix **[numBlocks_];
  for (int i=0; i<numBlocks_; i++) blockMatrix_[i] = new Epetra_CrsMatrix *[numBlocks_];
  for (int i=0; i<numBlocks_; i++)
    for (int j=0; j<numBlocks_; j++) blockMatrix_[i][j] = 0; // Initialize block operator pointers to zero

  blockGraph_ = new Epetra_Graph **[numBlocks_];
  for (int i=0; i<numBlocks_; i++) blockGraph_[i] = new Epetra_Graph *[numBlocks_];
  for (int i=0; i<numBlocks_; i++)
    for (int j=0; j<numBlocks_; j++) blockGraph_[i][j] = 0; // Initialize block graph pointers to zero

  blockLhs_ = new Epetra_Vector *[numBlocks_];
  for (int j=0; j<numBlocks_; j++) rowMap_[j] = 0; // Initialize Lhs pointers to zero

  blockRhs_ = new Epetra_Vector *[numBlocks_];
  for (int i=0; i<numBlocks_; i++) rowMap_[i] = 0; // Initialize Rhs pointers to zero

  rowMaps_ = new Epetra_Map *[numBlocks_];
  for (int i=0; i<numBlocks_; i++) rowMap_[i] = 0; // Initialize rowMap pointers to zero 

  colMaps_ = new Epetra_Map *[numBlocks_];
  for (int j=0; j<numBlocks_; j++) colMap_[j] = 0; // Initialize colMap pointers to zero

  return;
}

//=============================================================================
dft_SolverManager::~dft_SolverManager() {
 
  
  for (int i=0; i<numBlocks_; i++)
    for (int j=0; j<numBlocks_; j++) if (blockMatrix_[i][j]!=0) delete blockMatrix_[i][j];
  for (int i=0; i<numBlocks_; i++) delete [] blockMatrix_[i];
  delete [] blockMatrix_;

  for (int i=0; i<numBlocks_; i++)
    for (int j=0; j<numBlocks_; j++) if (blockGraph_[i][j]!=0) delete blockGraph_[i][j];
  for (int i=0; i<numBlocks_; i++) delete [] blockGraph_[i];
  delete [] blockGraph_;

  for (int j=0; j<numBlocks_; j++) delete blockLhs_[j];
  delete [] blockLhs_;

  for (int i=0; i<numBlocks_; i++) delete blockRhs_[i];
  delete [] blockRhs_;

  for (int i=0; i<numBlocks_; i++) delete rowMap_[i];
  delete [] rowMaps_;

  for (int j=0; j<numBlocks_; j++) delete colMap_[j];
  delete [] colMaps_;

  return;
}

//=============================================================================
int dft_SolverManager::finalizeBlockStructure() {
  
  if (isBlockStructureSet_) return(1); // Already been here, return warning
  // initialize graphs 
  for (int i=0; i<numBlocks_; i++) {
    blockLhs_ = new Epetra_Vector(rowMap_[i]);
    blockRhs_ = new Epetra_Vector(rowMap_[i]);
    for (int j=0; j<numBlocks_; j++) blockGraph_[i][j] = new Epetra_Graph(Copy, rowMap_[i], colMap_[j], 0);

    isBlockStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeGraphStructure() {
  
  if (isBlockStructureSet_) return(-1); // Block structure must be set
  if (isGraphStructureSet_) return(1); // Already been here, return warning

  // finalize graphs and optimize their storage, initialize matrix structure
  for (int i=0; i<numBlocks_; i++)
    for (int j=0; j<numBlocks_; j++) {
      EPETRA_CHK_ERR(blockGraph_[i][j]->FillComplete());
      EPETRA_CHK_ERR(blockGraph_[i][j]->OptimizeStorage());
      blockMatrix_[i][j] = new Epetra_CrsMatrix(Copy, blockGraph_[i][j]);
      EPETRA_CHK_ERR(blockMatrix_[i][j]->OptimizeStorage());
    }
  isGraphStructureSet_ = true;
  
  return(0);
}
//=============================================================================
int dft_SolverManager::initializeMatrixFill() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  
  // Zero out matrix Lhs and Rhs values
  for (int i=0; i<numBlocks_; i++) {
    EPETRA_CHK_ERR(blockLhs_->PutScalar(0.0));
    EPETRA_CHK_ERR(blockRhs_->PutScalar(0.0));
    for (int j=0; j<numBlocks_; j++) {
      EPETRA_CHK_ERR(blockMatrix_[i][j]->PutScalar(0.0));
    }
  }
  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeMatrixFill() {
  
  // Perform book-keeping after all values submitted
  for (int i=0; i<numBlocks_; i++)
    for (int j=0; j<numBlocks_; j++) {
      EPETRA_CHK_ERR(blockMatrix_[i][j]->FillComplete());
    }
  return(0);
}
//=============================================================================
int dft_SolverManager::setupSolver() {
  
  if (numBlocks_!=2) return(-1); // Can only do a two-block solver for now.
  // Perform book-keeping after all values submitted
  for (int i=0; i<numBlocks_; i++)
    for (int j=0; j<numBlocks_; j++) {
      EPETRA_CHK_ERR(blockMatrix_[i][j]->FillComplete());
    }
  return(0);
}
