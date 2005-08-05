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
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_TestForException.hpp"

//==============================================================================
dft_PolyA11_Epetra_Operator::dft_PolyA11_Epetra_Operator(const Epetra_Map & block1Map, const Epetra_Map * depNonLocalMap) 
  : block1Map_(block1Map),
    depNonLocalMap_(depNonLocalMap),
    matrix_(0),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true) {

  Label_ = "dft_HardSphereA11_Epetra_Operator";
  if (depNonLocalMap_!=0)
    matrix_ = new Epetra_CrsMatrix(Copy, *depNonLocalMap_, 0);
  

}
//==============================================================================
dft_PolyA11_Epetra_Operator::~dft_PolyA11_Epetra_Operator() {
  if (depNonLocalMap_!=0) delete matrix_;
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    if (matrix_!=0) 
      return(matrix_->PutScalar(0.0));
  
  return(0);
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {

  if (matrix_==0) return(0); // All ind NonLocal entries are diagonal values and 1, so we don't store them
  if (!depNonLocalMap_->MyGID(rowGID)) return(0);

  if (rowGID!=colGID) value = -value; // negate off-diagonal values to explicitly form inverse matrix

  if (firstTime_)
    return(matrix_->InsertGlobalValues(rowGID, 1, &value, &colGID));
  else
    return(matrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID));
  
  return(0);
}
//=============================================================================
int dft_PolyA11_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) 
    if (matrix_!=0) {
      matrix_->FillComplete(block1Map_, *depNonLocalMap_);
      matrix_->OptimizeStorage();
    }
  
  //  std::cout << *matrix_;
  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_PolyA11_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  Y=X; // We can safely do this

  if (matrix_!=0) { // Non-trivial depNonLocalMatrix

    int NumVectors = Y.NumVectors();
    int numMyElements = matrix_->NumMyRows();
    int offset = X.Length() - numMyElements; // We need to skip elements of Y associated with ind nonlocals
  

    double ** curY = new double *[NumVectors];
    double ** Yptr;
    
    Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
    for (int i=0; i<NumVectors; i++) curY[i] = Yptr[i]+offset;
    Epetra_MultiVector Ytmp(View, ownedMap_, curY, NumVectors); // Start Ytmp to view first numNodes elements of Y

    int ierr = matrix_->Multiply(false, X, Ytmp);
    delete [] curY;
  }
  return(ierr);
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
