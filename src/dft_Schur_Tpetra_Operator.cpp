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
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_Schur_Tpetra_Operator<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
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
  /* Used to capture matrix data
  LocalOrdinal nrows1 = A11->RowMap().NumGlobalElements();
  LocalOrdinal nrows2 = A22->RowMap().NumGlobalElements();
  LocalOrdinal ncols1 = A11->ColMap().NumGlobalElements();
  LocalOrdinal ncols2 = A22->ColMap().NumGlobalElements();
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
} //end constructor
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_Schur_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_Schur_Tpetra_Operator
()
{
} //end destructor
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_Schur_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{
  // Scalar normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in Schur Apply = " << normvalue << endl;
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
#endif

  // Apply (A22 - A21*inv(A11)*A12 to X

  RCP<MV > Y1 = rcp(new MV(A12_->getRangeMap(), X.getNumVectors()));
  RCP<MV > Y11 = rcp(new MV(A12_->getRangeMap(), X.getNumVectors()));
  RCP<MV > Y2 = rcp(new MV(A21_->getRangeMap(), X.getNumVectors()));

  Y.putScalar(STS::zero());
  A12op_->apply(X, *Y1);
  //Y1.NormInf(&normvalue);
  //cout << "Norm of Y1 in Schur Apply = " << normvalue << endl;
  //cout << *A12_ << endl;
  //exit(1);
  A11_->applyInverse(*Y1, *Y11);
  A21op_->apply(*Y11, *Y2);
  //Y2.NormInf(&normvalue);
  //cout << "Norm of Y2 in Schur Apply = " << normvalue << endl;
  A22_->apply(X, Y);
  Y.update(-STS::one(), *Y2, STS::one());

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in Schur Apply = " << normvalue << endl;
} //end Apply
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_Schur_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
ComputeRHS
(const MV& B1, const MV& B2, MV& B2S) const
{
  // Compute B2S =  B2 - A21*inv(A11)B1

  RCP<MV > Y1 = rcp(new MV(A11_->getDomainMap(), B1.getNumVectors()));
  A11_->applyInverse(B1, *Y1);
  A21op_->apply(*Y1, B2S);
  B2S.update(STS::one(), B2, -STS::one());
} //end ComputeRHS
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_Schur_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
ComputeX1
(const MV& B1, const MV& X2, MV& X1) const
{
  // Compute X1 =  inv(A11)(B1 - A12*X2)

  RCP<MV > Y1 = rcp(new MV(A12_->getRangeMap(), X1.getNumVectors()));
  A12op_->apply(X2, *Y1);
  Y1->update(STS::one(), B1, -STS::one());
  A11_->applyInverse(*Y1, X1);
} //end ComputeX1

//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_Schur_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
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

  // Apply (A22 - A21*inv(A11)*A12 to X

  RCP<MV > Y11 = rcp(new MV(A11_->getRangeMap(), X1.getNumVectors()));
  RCP<MV > Y12 = rcp(new MV(A12_->getRangeMap(), X1.getNumVectors()));
  RCP<MV > Y21 = rcp(new MV(A21_->getRangeMap(), X1.getNumVectors()));
  RCP<MV > Y22 = rcp(new MV(A22_->getRangeMap(), X1.getNumVectors()));

  Y1.putScalar(STS::zero());
  Y2.putScalar(STS::zero());
  A11_->apply(X1, *Y11);
  A12op_->apply(X2, *Y12);
  A21op_->apply(X1, *Y21);
  A22_->apply(X2, *Y22);
  Y1.update(STS::one(), *Y11, STS::one(), *Y12, STS::zero());
  Y2.update(STS::one(), *Y21, STS::one(), *Y22, STS::zero());
} //end ApplyGlobal
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_Schur_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
formSchurComplement
()
{
  if (A11invMatrix_.is_null() || A22Matrix_.is_null())
    {
      return;  // We cannot form S without the component matrices
    } //end if

  if (S_.get()==0)  // Form S
    {
      A11invA12_ = rcp(new MAT(A12_->getRowMap(), 0));
      A11invA12_->setObjectLabel("SchurComplement::A11invA12");
      A21A11invA12_ = rcp(new MAT(A21_->getRowMap(), 0));
      A21A11invA12_->setObjectLabel("SchurComplement::A21A11invA12");
      S_ = rcp(new MAT(A21_->getRowMap(), 0));
      S_->setObjectLabel("SchurComplement::S");
  } //end if

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
    S_->setAllToScalar(STS::zero());
  } //end if

  //Loop over rows and sum into
  size_t maxNumEntries = A22Matrix_->getGlobalMaxNumRowEntries();
  if(A21A11invA12_->getGlobalMaxNumRowEntries() > maxNumEntries)
    {
      maxNumEntries = A21A11invA12_->getGlobalMaxNumRowEntries();
    } //end if
  size_t NumEntries;
  Array<GlobalOrdinal> Indices(maxNumEntries);
  Array<MatScalar> Values(maxNumEntries);

  LocalOrdinal NumMyRows = S_->getNodeNumRows();
  LocalOrdinal Row, err;

  for( LocalOrdinal i = 0; i < NumMyRows; ++i )
  {
    Row = S_->getRowMap()->getGlobalElement(i);
    A22Matrix_->getGlobalRowCopy(Row, Indices, Values, NumEntries);
    if( sfilled ) //Sum In Values
    {
      S_->sumIntoGlobalValues( Row, Indices, Values);
    } //end if
    else
    {
      S_->insertGlobalValues( Row, Indices, Values);
    } //end else
    A21A11invA12_->getGlobalRowCopy( Row, Indices, Values, NumEntries);
    for( LocalOrdinal j = 0; j < NumEntries; ++j )
      {
	Values[j] = - Values[j];
      } //end for
    if( sfilled ) //Sum In Values
      {
	S_->sumIntoGlobalValues( Row, Indices, Values);
      } //end if
    else
      {
	S_->insertGlobalValues( Row, Indices, Values);
      } //end else
  } //end for

  if( !sfilled )
  {
    S_->fillComplete();
  } //end if

} //end getFormCompliment

TRAMONTO_INST_HELPER(dft_Schur_Tpetra_Operator)
