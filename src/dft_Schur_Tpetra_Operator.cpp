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

#include "dft_Schur_Tpetra_Operator.hpp"

//==============================================================================
template <class Scalar, class MatrixType>
dft_Schur_Tpetra_Operator<Scalar,MatrixType>::
dft_Schur_Tpetra_Operator
(RCP<APINV> A11, RCP<MAT> A12, RCP<MAT> A21, RCP<APINV> A22)
  : A11_(A11),
    A12_(A12),
    A21_(A21),
    A22_(A22),
    Label_(0)
{
  Label_ = "dft_Schur_Tpetra_Operator";
  A12op_ = rcp(new MMOP(A12_));
  A21op_ = rcp(new MMOP(A21_));
  A12rangeVec_ = rcp(new MV(A12_->getRangeMap(), 1));
  A12rangeVec2_ = rcp(new MV(A12_->getRangeMap(), 1));
  A21rangeVec_ = rcp(new MV(A21_->getRangeMap(), 1));
  A11rangeVec_ = rcp(new MV(A11_->getRangeMap(), 1));
  A11domainVec_ = rcp(new MV(A11_->getDomainMap(), 1));
  A22rangeVec_ = rcp(new MV(A22_->getRangeMap(), 1));

} //end constructor

//==============================================================================
template <class Scalar, class MatrixType>
dft_Schur_Tpetra_Operator<Scalar,MatrixType>::
~dft_Schur_Tpetra_Operator
()
{
} //end destructor

//==============================================================================
template <class Scalar, class MatrixType>
void
dft_Schur_Tpetra_Operator<Scalar,MatrixType>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
#endif

  Scalar ZERO = STS::zero();
  Scalar ONE = STS::one();

  // Implicitly apply the Schur complement (A22 - A21*inv(A11)*A12) to X
  Y.putScalar(ZERO);

  A12op_->apply(X, *A12rangeVec_);
  A11_->applyInverse(*A12rangeVec_, *A12rangeVec2_);
  A21op_->apply(*A12rangeVec2_, *A21rangeVec_);
  A22_->apply(X, Y);

  Y.update(-ONE, *A21rangeVec_, ONE);

} //end apply

//==============================================================================
template <class Scalar, class MatrixType>
void
dft_Schur_Tpetra_Operator<Scalar,MatrixType>::
ComputeRHS
(const MV& B1, const MV& B2, MV& B2S) const
{
  Scalar ONE = STS::one();

  // Compute B2S =  B2 - A21*inv(A11)B1
  A11_->applyInverse(B1, *A11domainVec_);
  A21_->apply(*A11domainVec_, B2S);
  B2S.update(ONE, B2, -ONE);
} //end ComputeRHS

//==============================================================================
template <class Scalar, class MatrixType>
void
dft_Schur_Tpetra_Operator<Scalar,MatrixType>::
ComputeX1
(const MV& B1, const MV& X2, MV& X1) const
{
  Scalar ONE = STS::one();

  // Compute X1 =  inv(A11)(B1 - A12*X2)
  A12op_->apply(X2, *A12rangeVec_);
  A12rangeVec_->update(ONE, B1, -ONE);
  A11_->applyInverse(*A12rangeVec_, X1);

} //end ComputeX1

//==============================================================================
template <class Scalar, class MatrixType>
void
dft_Schur_Tpetra_Operator<Scalar,MatrixType>::
ApplyGlobal
(const MV& X1, const MV& X2, MV& Y1, MV& Y2) const
{

  TEUCHOS_TEST_FOR_EXCEPT(Y1.getNumVectors()!=X1.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X1.getMap()->isSameAs(*A11_->getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!X2.getMap()->isSameAs(*A22_->getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y1.getMap()->isSameAs(*A11_->getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y2.getMap()->isSameAs(*A22_->getRangeMap()));
#endif

  Scalar ZERO = STS::zero();
  Scalar ONE = STS::one();

  // Apply (A22 - A21*inv(A11)*A12 to X
  Y1.putScalar(ZERO);
  Y2.putScalar(ZERO);

  A11_->apply(X1, *A11rangeVec_);
  A12op_->apply(X2, *A12rangeVec_);
  A21op_->apply(X1, *A21rangeVec_);
  A22_->apply(X2, *A22rangeVec_);
  Y1.update(ONE, *A11rangeVec_, ONE, *A12rangeVec_, ZERO);
  Y2.update(ONE, *A21rangeVec_, ONE, *A22rangeVec_, ZERO);
} //end ApplyGlobal

//==============================================================================
template <class Scalar, class MatrixType>
void
dft_Schur_Tpetra_Operator<Scalar,MatrixType>::
formSchurComplement
()
{
  if (A11invMatrix_.is_null() || A22Matrix_.is_null())
    {
      return;  // We cannot form S without the component matrices
    }

  if (S_.get()==0)  // Form S
    {
      A11invA12_ = rcp(new MAT(A12_->getRowMap(), 0));
      A11invA12_->setObjectLabel("SchurComplement::A11invA12");
      A21A11invA12_ = rcp(new MAT(A21_->getRowMap(), 0));
      A21A11invA12_->setObjectLabel("SchurComplement::A21A11invA12");
      S_ = rcp(new MAT(A21_->getRowMap(), 0));
      S_->setObjectLabel("SchurComplement::S");
  }

  bool A11invA12filled = A11invA12_->isFillComplete();
  bool A21A11invA12filled = A21A11invA12_->isFillComplete();

  if (A11invA12filled)
    A11invA12_->resumeFill();
  if (A21A11invA12filled)
    A11invA12_->resumeFill();

  // Compute inv(A11)*A12
    Tpetra::MatrixMatrix::Multiply(*A11invMatrix_, false, *A12_, false, *A11invA12_);
  // Compute A21A11invA12
    Tpetra::MatrixMatrix::Multiply(*A21_, false, *A11invA12_, false, *A21A11invA12_);

  // Finally compute S = A22 - A21A11invA12 manually

  //Initialize if S already filled
  bool sfilled = S_->isFillComplete();
  if (sfilled)
  {
    S_->resumeFill();
    S_->setAllToScalar(STMS::zero());
  }

  //Loop over rows and sum into
  size_t maxNumEntries = A22Matrix_->getGlobalMaxNumRowEntries();
  if(A21A11invA12_->getGlobalMaxNumRowEntries() > maxNumEntries)
    {
      maxNumEntries = A21A11invA12_->getGlobalMaxNumRowEntries();
    }
  size_t NumEntries;
  Array<GlobalOrdinal> Indices(maxNumEntries);
  Array<MatScalar> Values(maxNumEntries);

  LocalOrdinal NumMyRows = S_->getNodeNumRows();
  LocalOrdinal Row;

  for( LocalOrdinal i = OTLO::zero(); i < NumMyRows; ++i )
  {
    Row = S_->getRowMap()->getGlobalElement(i);
    A22Matrix_->getGlobalRowCopy(Row, Indices, Values, NumEntries);
    if( sfilled ) //Sum In Values
    {
      S_->sumIntoGlobalValues( Row, Indices, Values);
    }
    else
    {
      S_->insertGlobalValues( Row, Indices, Values);
    }
    A21A11invA12_->getGlobalRowCopy( Row, Indices, Values, NumEntries);
    for( LocalOrdinal j = OTLO::zero(); j < NumEntries; ++j )
    {
      Values[j] = - Values[j];
    }
    if( sfilled ) //Sum In Values
    {
      S_->sumIntoGlobalValues( Row, Indices, Values);
    }
    else
    {
      S_->insertGlobalValues( Row, Indices, Values);
    }
  }

  if( !sfilled )
  {
    S_->fillComplete();
  }

} //end getFormCompliment

TRAMONTO_INST_HELPER(dft_Schur_Tpetra_Operator)
