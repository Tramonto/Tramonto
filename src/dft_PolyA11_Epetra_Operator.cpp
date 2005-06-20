//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "dft_PolyA11_Epetra_Operator.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_InSerialDenseVector.h"

//==============================================================================
dft_PolyA11_Epetra_Operator::dft_PolyA11_Epetra_Operator(const Epetra_Map & ownedMap, int numBeads) 
  : ownedMap_(ownedMap),
    numBlocks_(2*numBeads),
    matrix_(0),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true) {

  Label_ = "dft_PolyA11_Epetra_Operator";
  matrix_ = Teuchos::rcp(new Teuchos::RefCountPtr<Epetra_CrsMatrix>[numBlocks_]);

  int numGlobalElements = ownedMap.NumGlobalElements();
  int numMyElements = ownedMap.NumMyElements();
  int * ownedGIDs = ownedMap.MyGlobalElements();
  Epetra_IntSerialDenseVector allGIDs(numBlocks_*numMyElements);
  int * ptr = allGIDs.Values();
  for (int i=0; i<numBlocks_; i++) {
    matrix_[i] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *ownedMap_.get(), 0));
    int offset = i*numGlobalElements;
    for (int j=0; j<numMyElements; j++) *ptr++ = offset+ownedGIDs[j];
  }
  vectorMap_ = Teuchos::rcp(new Epetra_Map(-1, numBlocks_*numMyElements, allGIDs.Values(), 0, ownedMap_.Comm()));
  
    
}
//==============================================================================
dft_PolyA11_Epetra_Operator::~dft_PolyA11_Epetra_Operator() {
}
//=============================================================================
int dft_SolverManager::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    for (int i=0; i<numBlocks_; i++)
      matrix_[i]->PutScalar(0.0);
  
  return(0);
}
//=============================================================================
int dft_SolverManager::insertMatrixValue(int ownedPhysicsID, int ownedNode, int colGID, double value) {

  if (firstTime_)
    matrix_[ownedPhysicsID]->InsertGlobalValues(ownedNode, 1, &value, &colGID);
  else
    matrix_[ownedPhysicsID]->SumIntoGlobalValues(ownedNode, 1, &value, &colGID);
  
  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  for (int i=0; i<numBlocks_; i++) {
    matrix_[i]->FillComplete(ownedMap_,vectorMap_.get());
    matrix_[i]->OptimizeStorage();
  }

  /*  for (int i=0; i<numBlocks_; i++) {
      std::cout << matrix_[i].get();
  */
  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  if (!X.Map().SameAs(OperatorDomainMap())) abort();  // These aborts should be handled as int return codes.
  if (!Y.Map().SameAs(OperatorRangeMap())) abort();
  if (Y.NumVectors()!=X.NumVectors()) abort();

  // Apply (A22 - A21*inv(A11)*A12 to X

  Epetra_MultiVector Y1(A12_->RangeMap(), X.NumVectors());
  Epetra_MultiVector Y11(A12_->RangeMap(), X.NumVectors());
  Epetra_MultiVector Y2(A21_->RangeMap(), X.NumVectors());
 
  A12_->Apply(X, Y1);
  A11_->Apply(Y1, Y11);
  A21_->Apply(Y11, Y2);
  A22_->Apply(X, Y);
  Y.Update(-1.0, Y2, 1.0);
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::ComputeRHS(const Epetra_MultiVector& B1, const Epetra_MultiVector& B2, 
					  Epetra_MultiVector& B2S) const {


  // Compute B2S =  B2 - A21*inv(A11)B1

  Epetra_MultiVector Y1(A11_->DomainMap(), B1.NumVectors());
 
  A11_->Apply(B1, Y1);
  A21_->Apply(Y1, B2S);
  B2S.Update(1.0, B2, -1.0);
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::ComputeX1(const Epetra_MultiVector& B1, const Epetra_MultiVector& X2, 
					 Epetra_MultiVector& X1) const {


  // Compute X1 =  inv(A11)(B1 - A12*X2)
 
  Epetra_MultiVector Y1(A12_->RangeMap(), X1.NumVectors());

  A12_->Apply(X2, Y1);
  Y1.Update(1.0, B1, -1.0);
  A11_->Apply(Y1, X1);
  return(0);
}
