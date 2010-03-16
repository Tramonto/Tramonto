//@HEADER
// ******************************************************************** 
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
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
  matrix_ = new Epetra_CrsMatrix*[numBlocks_-1];
  invDiagonal_ = new Epetra_Vector(block1Map);
  for (int i=0; i<numBlocks_-1; i++) {
    matrix_[i] = new Epetra_CrsMatrix(Copy, ownedMap, 0);
    matrix_[i]->SetLabel("PolyA11::matrix[i]");
    
  }
  return;
}
//==============================================================================
dft_PolyA11_Epetra_Operator::~dft_PolyA11_Epetra_Operator() {
  for (int i=0; i<numBlocks_-1; i++) if (matrix_[i]!=0) delete matrix_[i];
  delete [] matrix_;
  delete invDiagonal_;
  return;
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    for (int i=0; i<numBlocks_-1; i++)
      matrix_[i]->PutScalar(0.0);
    invDiagonal_->PutScalar(0.0);
  }
  return(0);
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::insertMatrixValue(int ownedPhysicsID, int ownedNode, int rowGID, int colGID, double value) {
  if (rowGID==colGID) {
    int locDiag = block1Map_.LID(colGID);
    (*invDiagonal_)[locDiag] += value;
    return(0);
  }
  else if (block1Map_.LID(colGID)> block1Map_.LID(rowGID)) {
    cout << "Encountered an illegal non-zero entry in dft_PolyA11_Epetra_Operator::insertMatrixValue." << endl
	 << "The A11 block cannot have nonzero terms in the upper diagonal." << endl
	 << "Input parameters:" << endl
	 << "  ownedPhysicsID = " << ownedPhysicsID << endl
	 << "  ownedNode      = " << ownedNode << endl
	 << "  rowGID         = " << rowGID << endl
	 << "  colGID         = " << colGID << endl
	 << "  block1Map_.LID(rowGID)         = " << block1Map_.LID(rowGID) << endl
	 << "  block1Map_.LID(colGID)         = " << block1Map_.LID(colGID) << endl
	 << "  value          = " << value << endl;
    return(-1);
  }

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
    matrix_[ownedPhysicsID-1]->SumIntoGlobalValues(ownedNode, 1, &value, &colGID);
  
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
  matrix_[curOwnedPhysicsID_-1]->InsertGlobalValues(curOwnedNode_, numEntries, values_.Values(), indices_.Values());

  curRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) 
    insertRow(); // Dump any remaining entries
    for (int i=0; i<numBlocks_-1; i++) {
      matrix_[i]->FillComplete(block1Map_, ownedMap_);
      matrix_[i]->OptimizeStorage();
      //cout << "PolyA11["<< i << "] Inf Norm = " << matrix_[i]->NormInf() << endl;
      //TEST_FOR_EXCEPT(!matrix_[i]->LowerTriangular());
    }
    invDiagonal_->Reciprocal(*invDiagonal_); // Invert diagonal values for faster ApplyInverse() method

    /*    for (int i=0; i<numBlocks_-1; i++) {
      std::cout << "matrix " << i << *matrix_[i];
      }*/
  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  //double normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in PolyA11 ApplyInverse = " << normvalue << endl;

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numMyElements = ownedMap_.NumMyElements();
  Epetra_MultiVector Ytmp(ownedMap_,NumVectors);

  Y=X; // We can safely do this

  double ** curYptr = new double *[NumVectors];
  double ** Yptr;
  Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
  for (int i=0; i<NumVectors; i++) curYptr[i] = Yptr[i];
  Epetra_MultiVector curY(View, ownedMap_, curYptr, NumVectors); // Start Ytmp to view first numNodes elements of Y

  double * curDiagptr;
  invDiagonal_->ExtractView(&curDiagptr); // Get pointer to first diagonal value
  Epetra_Vector curDiag(View, ownedMap_, curDiagptr); // View into first numNodes elements of diagonal

  curY.Multiply(1.0, curDiag, curY, 0.0); // Scale Y by the first block diagonal

  for (int i=0; i< numBlocks_-1; i++) { // Loop over block 1 through numBlocks (indexing 0 to numBlocks-1)
    // Update views of Y and diagonal blocks
    for (int j=0; j<NumVectors; j++) curYptr[j]+=numMyElements; // Increment pointers to next block
    curY.ResetView(curYptr); // Reset view to next block
    curDiagptr += numMyElements;
    curDiag.ResetView(curDiagptr); // Reset diagonal view to next block

    matrix_[i]->Multiply(false, Y, Ytmp); // Multiply block lower triangular block
    curY.Update(-1.0, Ytmp, 1.0); // curY = curX - Ytmp (Note that curX is in curY from initial copy Y = X)
    curY.Multiply(1.0, curDiag, curY, 0.0); // Scale Y by the first block diagonal
  }

  delete [] curYptr;

  // Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA11 ApplyInverse = " << normvalue << endl;
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  //double normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in PolyA11 Apply = " << normvalue << endl;

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numMyElements = ownedMap_.NumMyElements();

  double ** curYptr = new double *[NumVectors];
  double ** Yptr;

  Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
  for (int i=0; i<NumVectors; i++) curYptr[i] = Yptr[i];
  Epetra_MultiVector curY(View, ownedMap_, curYptr, NumVectors); // Start curY to view first numNodes elements of Y

  for (int i=0; i< numBlocks_-1; i++) {
    for (int j=0; j<NumVectors; j++) curYptr[j]+=numMyElements; // Increment pointers to next block
    curY.ResetView(curYptr); // Reset view to next block
    matrix_[i]->Multiply(false, X, curY); // This gives a result that is off-diagonal-matrix*X
  }

  Y.ReciprocalMultiply(1.0,*invDiagonal_, X, 1.0); // Add diagonal contribution

  delete [] curYptr;

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA11 Apply = " << normvalue << endl;

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
