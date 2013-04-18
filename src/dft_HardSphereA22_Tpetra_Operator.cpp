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

#include "dft_HardSphereA22_Tpetra_Operator.hpp"

//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_HardSphereA22_Tpetra_Operator
(const RCP<const MAP> & block2Map)
  : block2Map_(block2Map),
    densityOnDensityMatrix_(rcp(new VEC(block2Map))),
    densityOnDensityInverse_(rcp(new VEC(block2Map))),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true) {

  Label_ = "dft_HardSphereA22_Tpetra_Operator";
}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_HardSphereA22_Tpetra_Operator
()
{
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{
  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    densityOnDensityMatrix_->putScalar(STS::zero());
    densityOnDensityInverse_->putScalar(STS::zero());
  }

}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(GlobalOrdinal rowGID, GlobalOrdinal colGID, MatScalar value)
{
  TEUCHOS_TEST_FOR_EXCEPTION(rowGID!=colGID, std::logic_error, "Only diagonals are supposed to be entered.");

  // Storing this density block in a vector since it is diagonal
  densityOnDensityMatrix_->sumIntoLocalValue(block2Map_->getLocalElement(rowGID), value);

}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues()
{
  if (isLinearProblemSet_) {
    return; // nothing to do
  }

  // Form the inverse of densityOnDensityMatrix
  densityOnDensityInverse_->reciprocal(*densityOnDensityMatrix_);

  isLinearProblemSet_ = true;
  firstTime_ = false;

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const
{
  // Our algorithm is:
  // Y = D \ X

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
#endif

  Y.elementWiseMultiply(STS::one(), *densityOnDensityInverse_, X, STS::zero());

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
#endif

  Y.elementWiseMultiply(STS::one(), *densityOnDensityMatrix_, X, STS::one());

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
Check
(bool verbose) const
{

  RCP<VEC > x = rcp(new VEC(getDomainMap()));
  RCP<VEC > b = rcp(new VEC(getRangeMap()));
  x->randomize(); // Fill x with random numbers

  apply(*x, *b); // Forward operation
  applyInverse(*b, *b); // Reverse operation
  b->update(-STS::one(), *x, STS::one()); // Should be zero

  Scalar resid = b->norm2();

  if (verbose) {
    std::cout << "A22 self-check residual = " << resid << endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n");

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA22_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
formA22Matrix
()
{
  if (A22Matrix_ == Teuchos::null) {
    A22Matrix_ = rcp(new MAT(getRangeMap(), 1));
    A22Matrix_->setObjectLabel("HardSphereA22::A22Matrix");
  }
  LocalOrdinal numRows = getRangeMap()->getNodeNumElements();
  Teuchos::ArrayRCP<const Scalar> vectorValues = densityOnDensityMatrix_->get1dView();
  for (LocalOrdinal i=0; i<numRows; i++) {
    GlobalOrdinal row = A22Matrix_->getRowMap()->getGlobalElement(i);
    Scalar value = vectorValues[i];
    GlobalOrdinal col = row;
    Array<GlobalOrdinal> cols(1);
    Array<MatScalar> vals(1);
    cols[0] = col;
    vals[0] = Teuchos::as<MatScalar>(value);
    A22Matrix_->insertGlobalValues(row, cols, vals);
  }
  A22Matrix_->fillComplete();

}

TRAMONTO_INST_HELPER(dft_HardSphereA22_Tpetra_Operator)
