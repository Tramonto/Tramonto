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
    matrix_(0),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true) {

  Label_ = "dft_HardSphereA11_Epetra_Operator";
  if (depNonLocalMap_.NumGlobalElements()>0)
    matrix_ = new Epetra_CrsMatrix(Copy, depNonLocalMap_, 0);
    matrix_->SetLabel("HardSphere::A11::matrix");
  

}
//==============================================================================
dft_HardSphereA11_Epetra_Operator::~dft_HardSphereA11_Epetra_Operator() {
  if (matrix_!=0) delete matrix_;
}
//=============================================================================
int dft_HardSphereA11_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    if (matrix_!=0) 
      return(matrix_->PutScalar(0.0));
  
  return(0);
}
//=============================================================================
int dft_HardSphereA11_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {


 // All ind NonLocal entries are diagonal values and 1, so we don't store them
  if (matrix_==0) return(0); // No dependent nonlocal entries at all
  if (rowGID==colGID) return(0); // Don't keep diagonals
  if (!depNonLocalMap_.MyGID(rowGID)) return(0); // Isn't dependent nonlocal
  

  if (firstTime_)
    return(matrix_->InsertGlobalValues(rowGID, 1, &value, &colGID));
  else
    return(matrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID));
  
  return(0);
}
//=============================================================================
int dft_HardSphereA11_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) 
    if (matrix_!=0) {
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
  // | I  0 |
  // | B  I | 

  // We store only the X portion

  // The exact inverse is 
  // 
  // |  I  0 |
  // | -B  I | 


  if (matrix_==0) {
    Y=X;
    return(0);  // Nothing else to do
  }

  int NumVectors = Y.NumVectors();
  int numMyElements = matrix_->NumMyRows();
  int offset = Y.MyLength() - numMyElements; // We need to skip elements of Y associated with ind nonlocals
  assert(numMyElements==offset); // The two by two blocks should have the same dimension
  

  double ** curY = new double *[NumVectors];
  double ** Yptr;
  Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
  Epetra_MultiVector Y1(View, indNonLocalMap_, curY, NumVectors); // Start Y1 to view first block of Y elements
  for (int i=0; i<NumVectors; i++) curY[i] = Yptr[i]+offset;
  Epetra_MultiVector Y2(View, depNonLocalMap_, curY, NumVectors); // Start Y2 to view second block of Y elements
  
  double ** curX = new double *[NumVectors];
  double ** Xptr;
  X.ExtractView(&Xptr); // Get array of pointers to columns of X
  Epetra_MultiVector X1(View, indNonLocalMap_, curX, NumVectors); // Start X1 to view first block of X elements
  for (int i=0; i<NumVectors; i++) curX[i] = Xptr[i]+offset;
  Epetra_MultiVector X2(View, depNonLocalMap_, curX, NumVectors); // Start X2 to view second block of X elements
  
  int ierr = 0;
  if (&X[0]==&Y[0]) { // X and Y are the same
    // Y1 = X1 // Not needed
    Epetra_MultiVector Y2tmp(depNonLocalMap_, NumVectors); // Need space for multiplication result
    ierr = matrix_->Multiply(false, X1, Y2tmp); // Gives Y2tmp = B*X1
    if (inverse)
      Y2.Update(-1.0, Y2tmp, 1.0, X2, 0.0); // Gives us Y2 = X2 - B*X1
    else 
      Y2.Update(1.0, Y2tmp, 1.0, X2, 0.0); // Gives us Y2 = X2 + B*X1
  }
  else {
    Y1=X1;
    ierr = matrix_->Multiply(false, X1, Y2);// Gives Y2 = B*X1
    if (inverse) 
      Y2.Update(1.0, X2, -1.0); // Gives us Y2 = X2 - B*X1
    else
      Y2.Update(1.0, X2, 1.0); // Gives us Y2 = X2 + B*X1
  }
  delete [] curY;

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
