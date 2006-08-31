//@HEADER
// ******************************************************************** 
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER

#include "dft_HardSphereA22_Epetra_Operator.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Teuchos_TestForException.hpp"

//==============================================================================
dft_HardSphereA22_Epetra_Operator::dft_HardSphereA22_Epetra_Operator(const Epetra_Map & block2Map) 
  : block2Map_(block2Map),
    densityOnDensityMatrix_(Epetra_Vector(block2Map)),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true) {

  Label_ = "dft_HardSphereA22_Epetra_Operator";
}
//==============================================================================
dft_HardSphereA22_Epetra_Operator::~dft_HardSphereA22_Epetra_Operator() {
}
//=============================================================================
int dft_HardSphereA22_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    densityOnDensityMatrix_.PutScalar(0.0);
  }
  
  return(0);
}
//=============================================================================
int dft_HardSphereA22_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {
  
  assert(rowGID==colGID); // only diagonals are supposed to be entered

  densityOnDensityMatrix_[block2Map_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal

  return(0);
}
//=============================================================================
int dft_HardSphereA22_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do


  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_HardSphereA22_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  // Our algorithm is:
  // Y = D \ X


  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  Y.ReciprocalMultiply(1.0, densityOnDensityMatrix_, X, 0.0);

  return(0);
}
//==============================================================================
int dft_HardSphereA22_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  Y.Multiply(1.0, densityOnDensityMatrix_, X, 1.0);

  return(0);
}
//==============================================================================
int dft_HardSphereA22_Epetra_Operator::Check(bool verbose) const {

  Epetra_Vector x(OperatorDomainMap());
  Epetra_Vector b(OperatorRangeMap());
  x.Random(); // Fill x with random numbers
  Apply(x, b); // Forward operation
  ApplyInverse(b, b); // Reverse operation

  b.Update(-1.0, x, 1.0); // Should be zero

  double resid = 0.0;
  b.Norm2(&resid);

  if (verbose) 
    std::cout << "A22 self-check residual = " << resid << endl;

  if (resid > 1.0E-12) return(-1); // Bad residual
  return(0);
}
//==============================================================================
int dft_HardSphereA22_Epetra_Operator::formA22Matrix() {

  if (A22Matrix_.get()==0) {
    A22Matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, OperatorRangeMap(), 1)); 
    A22Matrix_->SetLabel("HardSphereA22::A22Matrix");

    int numRows = OperatorRangeMap().NumMyElements();
    for (int i=0; i<numRows; i++) {
      int row = A22Matrix_->GRID(i);
      double value = densityOnDensityMatrix_[i];
      int col = row;
      A22Matrix_->InsertGlobalValues(row, 1, &value, &col);
    }
  }
  else 
    A22Matrix_->ReplaceDiagonalValues(densityOnDensityMatrix_); 
  
  A22Matrix_->FillComplete();

  return(0);
}
