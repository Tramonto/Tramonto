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

#include "dft_BlockLinearProblem.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"

//==============================================================================
dft_BlockLinearProblem::dft_BlockLinearProblem(Epetra_Operator * A11inv, Epetra_Operator * A12, 
							     Epetra_Operator * A21, Epetra_Operator * A22) 
  : A11inv_(A11inv),
    A12_(A12),
    A21_(A21),
    A22_(A22),
    Label_(0) {

  Label_ = "dft_BlockLinearProblem";
}
//==============================================================================
dft_BlockLinearProblem::~dft_BlockLinearProblem() {
}
//==============================================================================
int dft_BlockLinearProblem::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  if (!X.Map().SameAs(OperatorDomainMap())){EPETRA_CHK_ERR(-1);}
  if (!Y.Map().SameAs(OperatorRangeMap())) {EPETRA_CHK_ERR(-2);}
  if (Y.NumVectors()!=X.NumVectors()) {EPETRA_CHK_ERR(-3);}

  // Apply (A22 - A21*inv(A11)*A12 to X

  Epetra_MultiVector Y1(A12_->RangeMap(), X.NumVectors());
  Epetra_MultiVector Y2(A21_->RangeMap(), X.NumVectors());
 
  A12_->Apply(X, Y1);
  A11inv_->Apply(Y1, Y1);
  A21_->Apply(Y1, Y2);
  A22_->Apply(X, Y);
  Y.Update(-1.0, Y2, 1.0);
  return(0);
}
//==============================================================================
double dft_BlockLinearProblem::ComputeGlobalResidual(const Epetra_Operator & A11, 
							    const Epetra_MultiVector& B1, const Epetra_MultiVector& B2, 
							    const Epetra_MultiVector& X1, const Epetra_MultiVector& X2) const {

  int NumVectors = X1.NumVectors();
  assert(NumVectors==X2.NumVectors());
  assert(NumVectors==B1.NumVectors());
  assert(NumVectors==B2.NumVectors());
  Epetra_MultiVector Y11(A11.RangeMap(), NumVectors);
  Epetra_MultiVector Y12(A12_->RangeMap(), NumVectors);
  Epetra_MultiVector Y21(A21_->RangeMap(), NumVectors);
  Epetra_MultiVector Y22(A22_->RangeMap(), NumVectors);

  EPETRA_CHK_ERR(A11.Apply(*X1, Y11));
  EPETRA_CHK_ERR(A12->Apply(*X2, Y12));
  EPETRA_CHK_ERR(A21->Apply(*X1, Y21));
  EPETRA_CHK_ERR(A22->Apply(*X2, Y22));
  
  // Compute Y11 = B1 - Y11 - Y12
  EPETRA_CHK_ERR(Y11.Update(1.0, B1, -1.0, Y12, -1.0));

  // Compute Y22 = B2 - Y22 - Y21
  EPETRA_CHK_ERR(Y22.Update(1.0, B2, -1.0, Y21, -1.0));

  Epetra_SerialDenseVector res1(NumVectors);
  Epetra_SerialDenseVector res2(NumVectors);

  EPETRA_CHK_ERR(Y11.NormInf(res1.Values())); // res1 and res2 now contain the inf-norms of each vector in the multivector
  EPETRA_CHK_ERR(Y22.NormInf(res2.Values()));

  double residual = EPETRA_MAX(res1.NormInf(), res2.NormInf()); // Now get the max over all max's.  This is the value we return

  return(residual);

//==============================================================================
int dft_BlockLinearProblem::ComputeGlobalRhs(const Epetra_Operator & A11, 
							    const Epetra_MultiVector& X1, const Epetra_MultiVector& X2, 
							    Epetra_MultiVector& B1, Epetra_MultiVector& B2) const {

  int NumVectors = X1.NumVectors();
  assert(NumVectors==X2.NumVectors());
  assert(NumVectors==B1.NumVectors());
  assert(NumVectors==B2.NumVectors());
  Epetra_MultiVector Y11(A11.RangeMap(), NumVectors);
  Epetra_MultiVector Y12(A12_->RangeMap(), NumVectors);
  Epetra_MultiVector Y21(A21_->RangeMap(), NumVectors);
  Epetra_MultiVector Y22(A22_->RangeMap(), NumVectors);

  EPETRA_CHK_ERR(A11.Apply(*X1, Y11));
  EPETRA_CHK_ERR(A12->Apply(*X2, Y12));
  EPETRA_CHK_ERR(A21->Apply(*X1, Y21));
  EPETRA_CHK_ERR(A22->Apply(*X2, Y22));
  
  // Compute B1 = Y11 + Y12
  EPETRA_CHK_ERR(B1.Update(1.0, Y11, 1.0, Y12, 0.0));

  // Compute B2 = Y22 + Y21
  EPETRA_CHK_ERR(B2.Update(1.0, Y22, 1.0, Y21, 0.0));

  return(0);
}
