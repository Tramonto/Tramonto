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

#include "dft_HardSphereA11_Tpetra_Operator.hpp"

//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_HardSphereA11_Tpetra_Operator
(const RCP<const MAP> & indNonLocalMap, const RCP<const MAP> & depNonLocalMap, const RCP<const MAP> & block1Map,
 RCP<ParameterList> parameterList)
  : indNonLocalMap_(indNonLocalMap),
    depNonLocalMap_(depNonLocalMap),
    block1Map_(block1Map),
    parameterList_(parameterList),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1) {

  Label_ = "dft_HardSphereA11_Tpetra_Operator";
  if (depNonLocalMap_->getGlobalNumElements()>0) {
    matrix_ = rcp(new MAT(depNonLocalMap_, 0));
    matrixOperator_ = rcp(new MMOP(matrix_));
    matrix_->setObjectLabel("HardSphere::A11::matrix");
  }


}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_HardSphereA11_Tpetra_Operator
()
{
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{
  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    if (matrix_!=Teuchos::null) {
      matrix_->resumeFill();
      matrix_->setAllToScalar(0.0);
    }
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal rowGID, LocalOrdinal colGID, MatScalar value) {

  Array<GlobalOrdinal> cols(1);
  cols[0] = colGID;
  Array<MatScalar> vals(1);
  vals[0] = value;

 // All ind NonLocal entries are diagonal values and 1, so we don't store them
  if (matrix_==Teuchos::null) {
    return; // No dependent nonlocal entries at all
  }
  if (rowGID==colGID) {
    return; // Don't keep diagonals
  }
  if (!depNonLocalMap_->isNodeGlobalElement(rowGID)) {
    return; // Isn't dependent nonlocal
  }

  if (firstTime_) {
    if (rowGID!=curRow_) {
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    }
    curRowValues_[colGID] += value;
  }
  else
    matrix_->sumIntoGlobalValues(rowGID, cols, vals);

}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  if (curRowValues_.empty())  {
    return;
  }
  size_t numEntries = curRowValues_.size();
  if (numEntries>indices_.size()) {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  }
  LocalOrdinal i=0;
  ITER pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  matrix_->insertGlobalValues(curRow_, indices_, values_);

  indices_.clear();
  values_.clear();
  curRowValues_.clear();

}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_) {
    return;
  }

  if (firstTime_) {
    insertRow(); // Dump any remaining entries
  }

  RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
  if (matrix_!=Teuchos::null) {
    matrix_->fillComplete(indNonLocalMap_, depNonLocalMap_, pl);
  }

  //  std::cout << *matrix_;
  isLinearProblemSet_ = true;
  firstTime_ = false;

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());

  // Matrix is of the form
  //
  // | -I  0 |
  // | B  -I |

  // We store only the X portion

  // The exact inverse is
  //
  // |  -I  0 |
  // | -B  -I |


  if (matrix_ == Teuchos::null) {
    Y.scale(-1.0, X); // Y = -X
    return;  // Nothing else to do
  }

  size_t NumVectors = Y.getNumVectors();
  size_t numMyElements = matrix_->getNodeNumRows();
  size_t offset = Y.getLocalLength() - numMyElements; // We need to skip elements of Y associated with ind nonlocals
  //Commented out 20-Nov-2006 assert(numMyElements==offset); // The two by two blocks should have the same dimension

  RCP<MV> Y1 = Y.offsetViewNonConst(indNonLocalMap_, 0);
  RCP<MV> Y2 = Y.offsetViewNonConst(depNonLocalMap_, offset);

  RCP<const MV> X1 = X.offsetView(indNonLocalMap_, 0);
  RCP<const MV> X2 = X.offsetView(depNonLocalMap_, offset);

  LocalOrdinal ierr = 0;
  if (&X.getVector(0)==&Y.getVector(0)) { // X and Y are the same
    RCP<MV> Y2tmp = rcp(new MV(depNonLocalMap_, NumVectors));
    matrixOperator_->apply(*X1, *Y2tmp);
    Y2->update(-1.0, *Y2tmp, -1.0, *X2, 0.0); // Gives us Y2 = -X2 - B*X1
    Y1->scale(-1.0, *X1);
  }
  else {
    Y1->scale(-1.0, *X1);
    matrixOperator_->apply(*X1, *Y2);
    Y2->update(-1.0, *X2, -1.0); // Gives us Y2 = -X2 - B*X1
  }

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());

  // Matrix is of the form
  //
  // | -I  0 |
  // | B  -I |

  // We store only the X portion

  if (matrix_ == Teuchos::null) {
    Y.scale(-1.0, X); // Y = -X
    return;  // Nothing else to do
  }

  size_t NumVectors = Y.getNumVectors();
  size_t numMyElements = matrix_->getNodeNumRows();
  size_t offset = Y.getLocalLength() - numMyElements; // We need to skip elements of Y associated with ind nonlocals
  //Commented out 20-Nov-2006 assert(numMyElements==offset); // The two by two blocks should have the same dimension

  RCP<MV> Y1 = Y.offsetViewNonConst(indNonLocalMap_, 0);
  RCP<MV> Y2 = Y.offsetViewNonConst(depNonLocalMap_, offset);

  RCP<const MV> X1 = X.offsetView(indNonLocalMap_, 0);
  RCP<const MV> X2 = X.offsetView(depNonLocalMap_, offset);

  LocalOrdinal ierr = 0;
  if (X.getVector(0)==Y.getVector(0)) { // X and Y are the same
    RCP<MV> Y2tmp = rcp(new MV(depNonLocalMap_, NumVectors));
    matrixOperator_->apply(*X1, *Y2tmp);
    Y2->update(1.0, *Y2tmp, -1.0, *X2, 0.0); // Gives us Y2 = -X2 + B*X1
    Y1->scale(-1.0, *X1);
  }
  else {
    Y1->scale(-1.0, *X1);
    matrixOperator_->apply(*X1, *Y2);
    Y2->update(-1.0, *X2, 1.0); // Gives us Y2 = -X2 + B*X1
  }

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
Check
(bool verbose) const
{
  RCP<VEC > x = rcp(new VEC(getDomainMap()));
  RCP<VEC > b = rcp(new VEC(getRangeMap()));

  x->randomize(); // Fill x with random numbers
  apply(*x, *b); // Forward operation
  applyInverse(*b, *b); // Reverse operation
  b->update(-1.0, *x, 1.0); // Should be zero

  Scalar absResid = b->norm2();
  Scalar normX = x->norm2();
  Scalar resid = absResid / normX;

  if (verbose) {
    std::cout << "A11 self-check residual = " << resid << std::endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n");

}
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
formA11invMatrix
()
{
  bool firstTime = false;
  if (A11invMatrix_ == Teuchos::null) {
    A11invMatrix_ = rcp(new MAT(getRangeMap(), 0));
    A11invMatrix_->setObjectLabel("HardSphereA11::A11invMatrix");
    firstTime = true;
  }
  else {
    A11invMatrix_->resumeFill();
    A11invMatrix_->setAllToScalar(0.0); // reset values
  }

  // insert -I for diagonal first
  LocalOrdinal numRows = getRangeMap()->getNodeNumElements();
  MatScalar value = -1.0;
  for (LocalOrdinal i=0; i<numRows; i++) {
    GlobalOrdinal row = A11invMatrix_->getRowMap()->getGlobalElement(i);
    GlobalOrdinal col = row;
    Array<GlobalOrdinal> cols(1);
    Array<MatScalar> vals(1);
    cols[0] = col;
    vals[0] = value;
    if (firstTime)
      A11invMatrix_->insertGlobalValues(row, cols, vals);
    else
      A11invMatrix_->sumIntoGlobalValues(row, cols, vals);
  }

  // Now insert lower triangle
  numRows = depNonLocalMap_->getNodeNumElements(); // number of rows in lower triangle
  if (numRows>0) {
    size_t numEntries;
    ArrayView<GlobalOrdinal> indices;
    ArrayView<MatScalar> values;
    for (LocalOrdinal i=0; i<numRows; i++) {
      GlobalOrdinal row = matrix_->getRowMap()->getGlobalElement(i);
      matrix_->getGlobalRowCopy( row, indices, values, numEntries );
      for( LocalOrdinal j = 0; j < numEntries; ++j )  {
	values[j] = - values[j];
      }
      if( firstTime ) {//Sum In Values
	A11invMatrix_->insertGlobalValues( row, indices, values );
      }
      else {
	A11invMatrix_->sumIntoGlobalValues( row, indices, values );
      }
    }

  }
  A11invMatrix_->fillComplete();

}
#if LINSOLVE_PREC == 0
// Use float
#if MIXED_PREC == 1
template class dft_HardSphereA11_Tpetra_Operator<float, float, int, int>;
#else
template class dft_HardSphereA11_Tpetra_Operator<float, float, int, int>;
#endif
#elif LINSOLVE_PREC == 1
// Use double
#if MIXED_PREC == 1
template class dft_HardSphereA11_Tpetra_Operator<double, float, int, int>;
#else
template class dft_HardSphereA11_Tpetra_Operator<double, double, int, int>;
#endif
#elif LINSOLVE_PREC == 2
// Use double double
#if MIXED_PREC == 1
template class dft_HardSphereA11_Tpetra_Operator<dd_real, double, int, int>;
#else
template class dft_HardSphereA11_Tpetra_Operator<dd_real, dd_real, int, int>;
#endif
#elif LINSOLVE_PREC == 3
// Use quad double
#if MIXED_PREC == 1
template class dft_HardSphereA11_Tpetra_Operator<qd_real, dd_real, int, int>;
#else
template class dft_HardSphereA11_Tpetra_Operator<qd_real, qd_real, int, int>;
#endif
#endif
