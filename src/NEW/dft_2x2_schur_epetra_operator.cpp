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

#include "dft_2x2_schur_epetra_operator.h"
#include "Epetra_Operator.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"

//==============================================================================
dft_2x2_schur_epetra_operator::dft_2x2_schur_epetra_operator(Epetra_Operator * A11inv, Epetra_Operator * A12, 
							     Epetra_Operator * A21, Epetra_Operator * A22) 
  : A11inv_(A11inv),
    A12_(A12),
    A21_(A21),
    A22_(A22),
    Label_(0) {

  Label_ = "dft_2x2_schur_epetra_operator";
}
//==============================================================================
dft_2x2_schur_epetra_operator::~dft_2x2_schur_epetra_operator() {
}
//==============================================================================
int dft_2x2_schur_epetra_operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


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
double dft_2x2_schur_epetra_operator::ComputeGlobalResidual(const Epetra_Operator & A11, 
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
int dft_2x2_schur_epetra_operator::ComputeGlobalRhs(const Epetra_Operator & A11, 
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

//==============================================================================
int dft_2x2_schur_epetra_operator::ComputeRHS(const Epetra_MultiVector& B1, const Epetra_MultiVector& B2, 
					  Epetra_MultiVector& B2S) const {


  // Compute B2S =  B2 - A21*inv(A11)B1

  Epetra_MultiVector Y1(A11inv_->DomainMap(), B1.NumVectors());
 
  A11inv_->Apply(B1, Y1);
  A21_->Apply(Y1, B2S);
  B2S.Update(1.0, B2, -1.0);
  return(0);
}
//==============================================================================
int dft_2x2_schur_epetra_operator::ComputeX1(const Epetra_MultiVector& B1, const Epetra_MultiVector& X2, 
					 Epetra_MultiVector& X1) const {


  // Compute X1 =  inv(A11)(B1 - A12*X2)
 
  A12_->Apply(X2, X1);
  X1.Update(1.0, B1, -1.0);
  A11inv_->Apply(X1, X1);
  return(0);
}
