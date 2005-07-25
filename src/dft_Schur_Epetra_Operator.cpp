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

#include "dft_Schur_Epetra_Operator.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Teuchos_TestForException.hpp"

//==============================================================================
dft_Schur_Epetra_Operator::dft_Schur_Epetra_Operator(Epetra_Operator * A11, Epetra_CrsMatrix * A12, 
						     Epetra_CrsMatrix * A21, Epetra_Operator * A22) 
  : A11_(A11),
    A12_(A12),
    A21_(A21),
    A22_(A22),
    Label_(0) {

  Label_ = "dft_Schur_Epetra_Operator";
  /* Used to capture matrix data
  int nrows1 = A11->RowMap().NumGlobalElements();
  int nrows2 = A22->RowMap().NumGlobalElements();
  int ncols1 = A11->ColMap().NumGlobalElements();
  int ncols2 = A22->ColMap().NumGlobalElements();
  Epetra_Map Row1Map(-1, nrows1, 0, A11->Comm());
  Epetra_Map Row2Map(-1, nrows2, nrows1, A22->Comm());
  Epetra_Map Col1Map(-1, ncols1, 0, A11->Comm());
  Epetra_Map Col2Map(-1, ncols2, ncols1, A22->Comm());

  assert(A11->ReplaceRowMap(Row1Map)==0);
  assert(A11->ReplaceColMap(Col1Map)==0);
  assert(A12->ReplaceRowMap(Row1Map)==0);
  assert(A12->ReplaceColMap(Col2Map)==0);
  assert(A21->ReplaceRowMap(Row2Map)==0);
  assert(A21->ReplaceColMap(Col1Map)==0);
  assert(A22->ReplaceRowMap(Row2Map)==0);
  assert(A22->ReplaceColMap(Col2Map)==0);

  EpetraExt::RowMatrixToMatrixMarketFile( "HardSphereA21.dat", *A21, "HardSphere_A21", 
					  "The 2,1 block of HardSphere Problem with Explicit non-local densities",
					  true);
  EpetraExt::RowMatrixToMatrixMarketFile( "HardSphereA12.dat", *A12, "HardSphere_A12", 
					  "The 1,2 block of HardSphere Problem with Explicit non-local densities",
					  true);
  EpetraExt::RowMatrixToMatrixMarketFile( "HardSphereA11.dat", *A11, "HardSphere_A11", 
					  "The 1,1 block of HardSphere Problem with Explicit non-local densities",
					  true);
  EpetraExt::RowMatrixToMatrixMarketFile( "HardSphereA22.dat", *A22, "HardSphere_A22", 
					  "The 2,2 block of HardSphere Problem with Explicit non-local densities",
					  true);
  abort();
  */
}
//==============================================================================
dft_Schur_Epetra_Operator::~dft_Schur_Epetra_Operator() {
}
//==============================================================================
int dft_Schur_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {


  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  // Apply (A22 - A21*inv(A11)*A12 to X

  Epetra_MultiVector Y1(A12_->RangeMap(), X.NumVectors());
  Epetra_MultiVector Y11(A12_->RangeMap(), X.NumVectors());
  Epetra_MultiVector Y2(A21_->RangeMap(), X.NumVectors());
 
  Y.PutScalar(0.0);
  A12_->Apply(X, Y1);
  A11_->ApplyInverse(Y1, Y11);
  A21_->Apply(Y11, Y2);
  A22_->Apply(X, Y);
  Y.Update(-1.0, Y2, 1.0);
  return(0);
}
//==============================================================================
int dft_Schur_Epetra_Operator::ComputeRHS(const Epetra_MultiVector& B1, const Epetra_MultiVector& B2, 
					  Epetra_MultiVector& B2S) const {


  // Compute B2S =  B2 - A21*inv(A11)B1

  Epetra_MultiVector Y1(A11_->OperatorDomainMap(), B1.NumVectors());
 
  A11_->ApplyInverse(B1, Y1);
  A21_->Apply(Y1, B2S);
  B2S.Update(1.0, B2, -1.0);
  return(0);
}
//==============================================================================
int dft_Schur_Epetra_Operator::ComputeX1(const Epetra_MultiVector& B1, const Epetra_MultiVector& X2, 
					 Epetra_MultiVector& X1) const {


  // Compute X1 =  inv(A11)(B1 - A12*X2)
 
  Epetra_MultiVector Y1(A12_->RangeMap(), X1.NumVectors());

  A12_->Apply(X2, Y1);
  Y1.Update(1.0, B1, -1.0);
  A11_->ApplyInverse(Y1, X1);
  return(0);
}
//==============================================================================
int dft_Schur_Epetra_Operator::ApplyGlobal(const Epetra_MultiVector& X1, const Epetra_MultiVector& X2, 
					   Epetra_MultiVector& Y1, Epetra_MultiVector& Y2) const {


  TEST_FOR_EXCEPT(!X1.Map().SameAs(A11_->OperatorDomainMap()));
  TEST_FOR_EXCEPT(!X2.Map().SameAs(A22_->OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y1.Map().SameAs(A11_->OperatorRangeMap()));
  TEST_FOR_EXCEPT(!Y2.Map().SameAs(A22_->OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y1.NumVectors()!=X1.NumVectors());

  // Apply (A22 - A21*inv(A11)*A12 to X

  Epetra_MultiVector Y11(A11_->OperatorRangeMap(), X1.NumVectors());
  Epetra_MultiVector Y12(A12_->OperatorRangeMap(), X1.NumVectors());
  Epetra_MultiVector Y21(A21_->OperatorRangeMap(), X1.NumVectors());
  Epetra_MultiVector Y22(A22_->OperatorRangeMap(), X1.NumVectors());
 
  Y1.PutScalar(0.0);
  Y2.PutScalar(0.0);
  A11_->Apply(X1, Y11);
  A12_->Apply(X2, Y12);
  A21_->Apply(X1, Y21);
  A22_->Apply(X2, Y22);
  Y1.Update(1.0, Y11, 1.0, Y12, 0.0);
  Y2.Update(1.0, Y21, 1.0, Y22, 0.0);
  return(0);
}
//==============================================================================
int dft_Schur_Epetra_Operator::CheckA11(bool verbose) const {

  Epetra_Vector x(A11_->OperatorDomainMap());
  Epetra_Vector b(A11_->OperatorRangeMap());
  x.Random(); // Fill x with random numbers
  A11_->Apply(x, b); // Forward operation
  A11_->ApplyInverse(b, b); // Reverse operation

  b.Update(-1.0, x, 1.0); // Should be zero

  double resid = 0.0;
  b.Norm2(&resid);

  if (verbose) 
    std::cout << "A11 self-check residual = " << resid << endl;

  if (resid > 1.0E-12) return(-1); // Bad residual
  return(0);
}
