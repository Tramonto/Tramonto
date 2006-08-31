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

#include "dft_schur_epetra_operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"

//==============================================================================
dft_schur_epetra_operator::dft_schur_epetra_operator(Epetra_CrsMatrix * A11, Epetra_CrsMatrix * A12, 
			    Epetra_CrsMatrix * A21, Epetra_CrsMatrix * A22) 
  : A11_(A11),
    A12_(A12),
    A21_(A21),
    A22_(A22),
    Label_(0) {

  Label_ = "dft_schur_epetra_operator";
}
//==============================================================================
dft_schur_epetra_operator::~dft_schur_epetra_operator() {
}
//==============================================================================
int dft_schur_epetra_operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  if (!X.Map().SameAs(OperatorDomainMap())) abort();  // These aborts should be handled as int return codes.
  if (!Y.Map().SameAs(OperatorRangeMap())) abort();
  if (Y.NumVectors()!=X.NumVectors()) abort();

  // Apply (A22 - A21*inv(A11)*A12 to X

  Epetra_MultiVector Y1(A12_->RangeMap(), X.NumVectors());
  Epetra_MultiVector Y2(A21_->RangeMap(), X.NumVectors());
 
  A12_->Apply(X, Y1);
  A11_->ApplyInverse(Y1, Y1);
  A21_->Apply(Y1, Y2);
  A22_->Apply(X, Y);
  Y.Update(-1.0, Y2, 1.0);
  return(0);
}
//==============================================================================
int dft_schur_epetra_operator::ComputeRHS(const Epetra_MultiVector& B1, const Epetra_MultiVector& B2, 
					  Epetra_MultiVector& B2S) const {


  // Compute B2S =  B2 - A21*inv(A11)B1

  Epetra_MultiVector Y1(A11_->DomainMap(), B1.NumVectors());
 
  A11_->ApplyInverse(B1, Y1);
  A21_->Apply(Y1, B2S);
  B2S.Update(1.0, B2, -1.0);
  return(0);
}
//==============================================================================
int dft_schur_epetra_operator::ComputeX1(const Epetra_MultiVector& B1, const Epetra_MultiVector& X2, 
					 Epetra_MultiVector& X1) const {


  // Compute X1 =  inv(A11)(B1 - A12*X2)
 
  A12_->Apply(X2, X1);
  X1.Update(1.0, B1, -1.0);
  A11_->ApplyInverse(X1, X1);
  return(0);
}
