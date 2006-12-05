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

#include "dft_PolyA11_Epetra_Operator.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_TestForException.hpp"

//==============================================================================
dft_PolyA11_Epetra_Operator::dft_PolyA11_Epetra_Operator(const Epetra_Map & ownedMap, const Epetra_Map & block1Map) 
  : ownedMap_(ownedMap),
    block1Map_(block1Map),
    numBlocks_(block1Map.NumMyElements()/ownedMap.NumMyElements()),
    matrix_(0),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1),
    curOwnedPhysicsID_(-1),
    curOwnedNode_(-1) {

  Label_ = "dft_PolyA11_Epetra_Operator";
  matrix_ = new Epetra_CrsMatrix*[numBlocks_];
  for (int i=0; i<numBlocks_; i++) {
    matrix_[i] = new Epetra_CrsMatrix(Copy, ownedMap, 0);
    matrix_[i]->SetLabel("PolyA11::matrix[i]");
  }

}
//==============================================================================
dft_PolyA11_Epetra_Operator::~dft_PolyA11_Epetra_Operator() {
  for (int i=0; i<numBlocks_; i++) if (matrix_[i]!=0) delete matrix_[i];
  delete [] matrix_;
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    for (int i=0; i<numBlocks_; i++)
      matrix_[i]->PutScalar(0.0);
  
  return(0);
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::insertMatrixValue(int ownedPhysicsID, int ownedNode, int rowGID, int colGID, double value) {

  if (rowGID!=colGID) value = -value; // negate off-diagonal values to simplify kernel calls

  if (firstTime_) {
    if (rowGID!=curRow_) { 
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
      curOwnedPhysicsID_ = ownedPhysicsID;
      curOwnedNode_ = ownedNode;
    }
    curRowValues_[colGID] += value;
  }
  else
    matrix_[ownedPhysicsID]->SumIntoGlobalValues(ownedNode, 1, &value, &colGID);
  
  return(0);
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::insertRow() {
  if (curRowValues_.empty()) return(0);
  int numEntries = curRowValues_.size();
  if (numEntries>indices_.Length()) {
    indices_.Resize(numEntries);
    values_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  matrix_[curOwnedPhysicsID_]->InsertGlobalValues(curOwnedNode_, numEntries, values_.Values(), indices_.Values());

  curRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) 
    insertRow(); // Dump any remaining entries
    for (int i=0; i<numBlocks_; i++) {
      matrix_[i]->FillComplete(block1Map_, ownedMap_);
      matrix_[i]->OptimizeStorage();
      //TEST_FOR_EXCEPT(!matrix_[i]->LowerTriangular());
    }
  
  /*  for (int i=0; i<numBlocks_; i++) {
      std::cout << *matrix_[i];
  */
  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numMyElements = ownedMap_.NumMyElements();

  Y=X; // We can safely do this

  double ** curY = new double *[NumVectors];
  double ** Yptr;

  Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
  for (int i=0; i<NumVectors; i++) curY[i] = Yptr[i];
  Epetra_MultiVector Ytmp(View, ownedMap_, curY, NumVectors); // Start Ytmp to view first numNodes elements of Y

  for (int i=0; i< numBlocks_; i++) {
    matrix_[i]->Multiply(false, Y, Ytmp);
    for (int j=0; j<NumVectors; j++) curY[j]+=numMyElements; // Increment pointers to next block
    Ytmp.ResetView(curY); // Reset view to next block
  }
  delete [] curY;
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numMyElements = ownedMap_.NumMyElements();

  double ** curY = new double *[NumVectors];
  double ** curX = new double *[NumVectors];
  double ** Yptr;
  double ** Xptr;

  Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
  for (int i=0; i<NumVectors; i++) curY[i] = Yptr[i];
  Epetra_MultiVector Ytmp(View, ownedMap_, curY, NumVectors); // Start Ytmp to view first numNodes elements of Y
  X.ExtractView(&Xptr); // Get array of pointers to columns of X
  for (int i=0; i<NumVectors; i++) curX[i] = Xptr[i];
  Epetra_MultiVector Xtmp(View, ownedMap_, curX, NumVectors); // Start Xtmp to view first numNodes elements of X

  for (int i=0; i< numBlocks_; i++) {
    matrix_[i]->Multiply(false, X, Ytmp); // This gives a result that is X - off-diagonal-matrix*X
    Ytmp.Update(-2.0, Xtmp, 1.0); // This gives a result of -X - off-diagonal-matrix*X
    Ytmp.Scale(-1.0); // Finally negate to get the desired result
    for (int j=0; j<NumVectors; j++) {
      curY[j]+=numMyElements; // Increment pointers to next block
      curX[j]+=numMyElements; // Increment pointers to next block
    }
    Ytmp.ResetView(curY); // Reset view to next block
    Xtmp.ResetView(curX); // Reset view to next block
  }
  delete [] curY;
  delete [] curX;

  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::Check(bool verbose) const {

  Epetra_Vector x(OperatorDomainMap());
  Epetra_Vector b(OperatorRangeMap());
  x.Random(); // Fill x with random numbers
  Apply(x, b); // Forward operation
  ApplyInverse(b, b); // Reverse operation

  b.Update(-1.0, x, 1.0); // Should be zero

  double resid = 0.0;
  b.Norm2(&resid);

  if (verbose) 
    std::cout << "A11 self-check residual = " << resid << endl;

  if (resid > 1.0E-12) return(-1); // Bad residual
  return(0);
}
