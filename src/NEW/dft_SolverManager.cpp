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
  : comm_(comm),
    blockMatrix_(0),
    blockLhs_(0),
    blockRhs_(0),
    rowMaps_(0),
    colMaps_(0),
    isStructureSet_(false) {
  // Allocate blockMatrix, blockGraph, blockLhs, blockRhs and map arrays
  blockMatrix_ = new Epetra_CrsMatrix **[numBlocks];
  for (int i=0; i<numBlocks; i++) blockMatrix_[i] = new Epetra_CrsMatrix *[numBlocks];
  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) blockMatrix_[i][j] = 0; // Initialize block operator pointers to zero

  blockGraph_ = new Epetra_Graph **[numBlocks];
  for (int i=0; i<numBlocks; i++) blockGraph_[i] = new Epetra_Graph *[numBlocks];
  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) blockGraph_[i][j] = 0; // Initialize block graph pointers to zero

  blockLhs_ = new Epetra_Vector *[numBlocks];
  for (int j=0; j<numBlocks; j++) rowMap_[j] = 0; // Initialize Lhs pointers to zero

  blockRhs_ = new Epetra_Vector *[numBlocks];
  for (int i=0; i<numBlocks; i++) rowMap_[i] = 0; // Initialize Rhs pointers to zero

  rowMaps_ = new Epetra_Map *[numBlocks];
  for (int i=0; i<numBlocks; i++) rowMap_[i] = 0; // Initialize rowMap pointers to zero 

  colMaps_ = new Epetra_Map *[numBlocks];
  for (int j=0; j<numBlocks; j++) colMap_[j] = 0; // Initialize colMap pointers to zero

  return;
}

//=============================================================================
dft_SolverManager::~dft_SolverManager() {
 
  
  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) delete blockMatrix_[i][j];
  for (int i=0; i<numBlocks; i++) delete [] blockMatrix_[i];
  delete [] blockMatrix_;

  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) delete blockGraph_[i][j];
  for (int i=0; i<numBlocks; i++) delete [] blockGraph_[i];
  delete [] blockGraph_;

  for (int j=0; j<numBlocks; j++) delete rowMap_[j];
  delete [] blockLhs_;

  for (int i=0; i<numBlocks; i++) delete rowMap_[i];
  delete [] blockRhs_;

  for (int i=0; i<numBlocks; i++) delete rowMap_[i];
  delete [] rowMaps_;

  for (int j=0; j<numBlocks; j++) delete colMap_[j];
  delete [] colMaps_;

  return;
}

//=============================================================================
int dft_SolverManager::finalizeBlockStructure() {
  
  // initialize graphs 
  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) blockGraph_[i][j] = new Epetra_Graph(Copy, rowMap_[i], colMap_[j], 0);
  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeGraphStructure() {
  
  // finalize graphs and optimize their storage, initialize matrix structure
  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) {
      EPETRA_CHK_ERR(blockGraph_[i][j]->FillComplete());
      EPETRA_CHK_ERR(blockGraph_[i][j]->OptimizeStorage());
      blockMatrix_[i][j] = new Epetra_CrsMatrix(Copy, blockGraph_[i][j]);
      EPETRA_CHK_ERR(blockMatrix_[i][j]->OptimizeStorage());
    }
  return(0);
}
//=============================================================================
int dft_SolverManager::initializeMatrixFill() {
  
  // Zero out matrix values
  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) {
      EPETRA_CHK_ERR(blockMatrix_[i][j]->PutScalar(0.0));
    }
  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeMatrixFill() {
  
  // Perform book-keeping after all values submitted
  for (int i=0; i<numBlocks; i++; i++)
    for (int j=0; j<numBlocks; j++) {
      EPETRA_CHK_ERR(blockMatrix_[i][j]->FillComplete());
    }
  return(0);
}
