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

#include "dft_HardSphereA11_Epetra_Operator.hpp"
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
dft_HardSphereA11_Epetra_Operator::dft_HardSphereA11_Epetra_Operator(const Epetra_Map & indNonLocalMap, const Epetra_Map & depNonLocalMap, const Epetra_Map & block1Map) 
  : indNonLocalMap_(indNonLocalMap),
    depNonLocalMap_(depNonLocalMap),
    block1Map_(block1Map),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1) {

  Label_ = "dft_HardSphereA11_Epetra_Operator";
  if (depNonLocalMap_.NumGlobalElements()>0) {
    matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, depNonLocalMap_, 0));
    matrix_->SetLabel("HardSphere::A11::matrix");
  }
  

}
//==============================================================================
dft_HardSphereA11_Epetra_Operator::~dft_HardSphereA11_Epetra_Operator() {
}
//=============================================================================
int dft_HardSphereA11_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    if (matrix_.get()!=0) 
      return(matrix_->PutScalar(0.0));
  
  return(0);
}
//=============================================================================
int dft_HardSphereA11_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {


 // All ind NonLocal entries are diagonal values and 1, so we don't store them
  if (matrix_.get()==0) return(0); // No dependent nonlocal entries at all
  if (rowGID==colGID) return(0); // Don't keep diagonals
  if (!depNonLocalMap_.MyGID(rowGID)) return(0); // Isn't dependent nonlocal
  

  if (firstTime_) {
    if (rowGID!=curRow_) { 
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    }
    curRowValues_[colGID] += value;
  }
  else
    return(matrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID));
  
  return(0);
}
//=============================================================================
int dft_HardSphereA11_Epetra_Operator::insertRow() {
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
  matrix_->InsertGlobalValues(curRow_, numEntries, values_.Values(), indices_.Values());

  curRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_HardSphereA11_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) 
    if (matrix_.get()!=0) {
      insertRow(); // Dump any remaining entries
      matrix_->FillComplete(indNonLocalMap_, depNonLocalMap_);
      matrix_->OptimizeStorage();
    }
  
  //  std::cout << *matrix_;
  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_HardSphereA11_Epetra_Operator::doApply(const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool inverse) const {


  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  // Matrix is of the form 
  // 
  // | -I  0 |
  // | B  -I | 

  // We store only the X portion

  // The exact inverse is 
  // 
  // |  -I  0 |
  // | -B  -I | 


  if (matrix_.get()==0) {
    Y.Scale(-1.0, X); // Y = -X
    return(0);  // Nothing else to do
  }

  int NumVectors = Y.NumVectors();
  int numMyElements = matrix_->NumMyRows();
  int offset = Y.MyLength() - numMyElements; // We need to skip elements of Y associated with ind nonlocals
  //Commented out 20-Nov-2006 assert(numMyElements==offset); // The two by two blocks should have the same dimension
  

  double ** curY = new double *[NumVectors];
  double ** Yptr;
  Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
  Epetra_MultiVector Y1(View, indNonLocalMap_, Yptr, NumVectors); // Start Y1 to view first block of Y elements
  for (int i=0; i<NumVectors; i++) curY[i] = Yptr[i]+offset;
  Epetra_MultiVector Y2(View, depNonLocalMap_, curY, NumVectors); // Start Y2 to view second block of Y elements
  
  double ** curX = new double *[NumVectors];
  double ** Xptr;
  X.ExtractView(&Xptr); // Get array of pointers to columns of X
  Epetra_MultiVector X1(View, indNonLocalMap_, Xptr, NumVectors); // Start X1 to view first block of X elements
  for (int i=0; i<NumVectors; i++) curX[i] = Xptr[i]+offset;
  Epetra_MultiVector X2(View, depNonLocalMap_, curX, NumVectors); // Start X2 to view second block of X elements
  
  int ierr = 0;
  if (&X[0]==&Y[0]) { // X and Y are the same
    Epetra_MultiVector Y2tmp(depNonLocalMap_, NumVectors); // Need space for multiplication result
    ierr = matrix_->Multiply(false, X1, Y2tmp); // Gives Y2tmp = B*X1
    if (inverse)
      Y2.Update(-1.0, Y2tmp, -1.0, X2, 0.0); // Gives us Y2 = -X2 - B*X1
    else 
      Y2.Update(1.0, Y2tmp, -1.0, X2, 0.0); // Gives us Y2 = -X2 + B*X1
    Y1.Scale(-1.0,X1); // Y1 = -X1 
  }
  else {
    Y1.Scale(-1.0, X1);
    ierr = matrix_->Multiply(false, X1, Y2);// Gives Y2 = B*X1
    if (inverse) 
      Y2.Update(-1.0, X2, -1.0); // Gives us Y2 = -X2 - B*X1
    else
      Y2.Update(-1.0, X2, 1.0); // Gives us Y2 = -X2 + B*X1
  }
  delete [] curY;
  delete [] curX;

  return(ierr);
}
//==============================================================================
int dft_HardSphereA11_Epetra_Operator::Check(bool verbose) const {

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
//==============================================================================
int dft_HardSphereA11_Epetra_Operator::formA11invMatrix() {

  bool firstTime = false;
  if (A11invMatrix_.get()==0) {
    A11invMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, OperatorRangeMap(), 0)); 
    A11invMatrix_->SetLabel("HardSphereA11::A11invMatrix");
    firstTime = true;
  }
  else A11invMatrix_->PutScalar(0.0); // reset values

  // insert -I for diagonal first

  int numRows = OperatorRangeMap().NumMyElements();
  double value = -1.0;
  for (int i=0; i<numRows; i++) {
    int row = A11invMatrix_->GRID(i);
    int col = row;
    if (firstTime)
      A11invMatrix_->InsertGlobalValues(row, 1, &value, &col);
    else
      A11invMatrix_->SumIntoGlobalValues(row, 1, &value, &col);
  }
  // Now insert lower triangle
  numRows = depNonLocalMap_.NumMyElements(); // number of rows in lower triangle
  if (numRows>0) {
    int MaxNumEntries = matrix_->MaxNumEntries();
    int NumEntries;
    int * Indices = new int[MaxNumEntries];
    double * Values = new double[MaxNumEntries];
    for (int i=0; i<numRows; i++) {
      int Row = matrix_->GRID(i);
      EPETRA_CHK_ERR( matrix_->ExtractGlobalRowCopy( Row, MaxNumEntries, NumEntries, Values, Indices ) );
      for( int j = 0; j < NumEntries; ++j ) Values[j] = - Values[j];
      if( firstTime ) {//Sum In Values
	int err = A11invMatrix_->InsertGlobalValues( Row, NumEntries, Values, Indices ); assert( err == 0 );
      }
      else {
	int err = A11invMatrix_->SumIntoGlobalValues( Row, NumEntries, Values, Indices ); assert( err == 0 || err == 1 );
      }
    }
    
    delete [] Indices;
    delete [] Values;
  }
  A11invMatrix_->FillComplete();
  return(0);
}
