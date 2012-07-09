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

#include "dft_PolyA11_Tpetra_Operator.hpp"

///==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA11_Tpetra_Operator
(const RCP<const MAP> & ownedMap, const RCP<const MAP> & block1Map)
  : ownedMap_(ownedMap),
    block1Map_(block1Map),
    numBlocks_(block1Map->getNodeNumElements()/ownedMap->getNodeNumElements()),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1),
    curOwnedPhysicsID_(-1),
    curOwnedNode_(-1)
{

  Label_ = "dft_PolyA11_Tpetra_Operator";

  invDiagonal_ = rcp(new VEC(block1Map));
  matrix_.resize(numBlocks_-1);
  for (LocalOrdinal i=0; i<numBlocks_-1; i++)
  {
    matrix_[i] = rcp(new MAT(ownedMap, 0));
    matrix_[i]->setObjectLabel("PolyA11::matrix[i]");
  } //end for
  return;
} //end constructor
//=============================================================================
// RN_20100326: A constructor that in addition takes a parameter list. This
// is for Krylov solvers used in the applyInverse method.
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA11_Tpetra_Operator
(const RCP<const MAP > & ownedMap, const RCP<const MAP > & block1Map,
 RCP<ParameterList> parameterList)
  : ownedMap_(ownedMap),
    block1Map_(block1Map),
    parameterList_(parameterList),
    numBlocks_(block1Map->getNodeNumElements()/ownedMap->getNodeNumElements()),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1),
    curOwnedPhysicsID_(-1),
    curOwnedNode_(-1)
{
  Label_ = "dft_PolyA11_Tpetra_Operator";

  invDiagonal_ = rcp(new VEC(block1Map));

  matrix_.resize(numBlocks_-1);
  for (LocalOrdinal i=0; i<numBlocks_-1; i++)
  {
    matrix_[i] = rcp(new MAT(ownedMap, 0));
    matrix_[i]->setObjectLabel("PolyA11::matrix[i]");
  } //end for

  return;
} //end constructor
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_PolyA11_Tpetra_Operator
()
{
  return;
} //end destructor
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{
  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
  {
    for (LocalOrdinal i=0; i<numBlocks_-1; i++)
    {
      matrix_[i]->resumeFill();
      matrix_[i]->setAllToScalar(0.0);
    } //end for

    invDiagonal_->putScalar(0.0);
  } //end if

} //end initializeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value)
{

  Array<GlobalOrdinal> cols(1);
  Array<Scalar> vals(1);
  cols[0] = colGID;
  vals[0] = value;
  if (rowGID==colGID)
  {
    LocalOrdinal locDiag = block1Map_->getLocalElement(colGID);
    invDiagonal_->sumIntoLocalValue(locDiag, value);
    return;
  } //end if
  TEUCHOS_TEST_FOR_EXCEPTION(block1Map_->getLocalElement(colGID)> block1Map_->getLocalElement(rowGID), std::runtime_error,
    std::cout << "Encountered an illegal non-zero entry in dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::insertMatrixValue." << std::endl
	 << "The A11 block cannot have nonzero terms in the upper diagonal." << std::endl
	 << "Input parameters:" << std::endl
	 << "  ownedPhysicsID = " << ownedPhysicsID << std::endl
	 << "  ownedNode      = " << ownedNode << std::endl
	 << "  rowGID         = " << rowGID << std::endl
	 << "  colGID         = " << colGID << std::endl
	 << "  block1Map_.LID(rowGID)         = " << block1Map_->getLocalElement(rowGID) << std::endl
	 << "  block1Map_.LID(colGID)         = " << block1Map_->getLocalElement(colGID) << std::endl
	 << "  value          = " << value << std::endl);

  if (firstTime_)
  {
    if (rowGID!=curRow_)
    {
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
      curOwnedPhysicsID_ = ownedPhysicsID;
      curOwnedNode_ = ownedNode;
    } //end if
    curRowValues_[colGID] += value;
  } //end if
  else
  {
    matrix_[ownedPhysicsID-1]->sumIntoGlobalValues(ownedNode, cols, vals);

  } //end else
} //end insertMatrixValue
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  if (curRowValues_.empty())
  {
    return;
  } //end if
  size_t numEntries = curRowValues_.size();
  if (numEntries>indices_.size())
  {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  } //end if
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos)
  {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  } //end for

  matrix_[curOwnedPhysicsID_-1]->insertGlobalValues(curOwnedNode_, indices_, values_);

  indices_.clear();
  values_.clear();
  curRowValues_.clear();
} //end insertRow
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
  {
    return; // nothing to do
  } //end if

  if (firstTime_)
  {
    insertRow(); // Dump any remaining entries
  } //end if
  for (LocalOrdinal i=0; i<numBlocks_-1; i++)
  {
    matrix_[i]->fillComplete(block1Map_, ownedMap_);
    //cout << "PolyA11["<< i << "] Inf Norm = " << matrix_[i]->NormInf() << endl;
    //TEUCHOS_TEST_FOR_EXCEPT(!matrix_[i]->LowerTriangular());
  } //end for
  invDiagonal_->reciprocal(*invDiagonal_); // Invert diagonal values for faster applyInverse() method

  /*
  for (LocalOrdinal i=0; i<numBlocks_-1; i++)
  {
    std::cout << "matrix " << i << *matrix_[i];
  } //end for
  */
  isLinearProblemSet_ = true;
  firstTime_ = false;
} //end finalizeProblemValues
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  printf("\n\n\n\ndft_PolyA11_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  size_t NumVectors = Y.getNumVectors();
  size_t numMyElements = ownedMap_->getNodeNumElements();
  RCP<MV > Ytmp = rcp(new MV(ownedMap_,NumVectors));

  Y=X; // We can safely do this

  RCP<MV > curY = Y.offsetViewNonConst(ownedMap_, 0);
  // Start Ytmp to view first numNodes elements of Y

  RCP<VEC> diagVec = invDiagonal_->offsetViewNonConst(ownedMap_, 0)->getVectorNonConst(0);

  curY->elementWiseMultiply(1.0, *diagVec, *curY, 0.0); // Scale Y by the first block diagonal

  // Loop over block 1 through numBlocks (indexing 0 to numBlocks-1)
  for (LocalOrdinal i=0; i< numBlocks_-1; i++)
  {
    // Update views of Y and diagonal blocks
    //for (LocalOrdinal j=0; j<NumVectors; j++)
    curY = Y.offsetViewNonConst(ownedMap_, (i+1)*numMyElements);

    diagVec = invDiagonal_->offsetViewNonConst(ownedMap_, (i+1)*numMyElements)->getVectorNonConst(0);

    matrix_[i]->apply(Y, *Ytmp); // Multiply block lower triangular block
    curY->update(-1.0, *Ytmp, 1.0); // curY = curX - Ytmp (Note that curX is in curY from initial copy Y = X)
    curY->elementWiseMultiply(1.0, *diagVec, *curY, 0.0); // Scale Y by the first block diagonal
  } //end for
} //end applyInverse
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
  size_t NumVectors = Y.getNumVectors();
  size_t numMyElements = ownedMap_->getNodeNumElements();

  RCP<MV > curY = Y.offsetViewNonConst(ownedMap_, 0);
  // Start curY to view first numNodes elements of Y

  for (LocalOrdinal i=0; i< numBlocks_-1; i++) {
    curY = Y.offsetViewNonConst(ownedMap_, (i+1)*numMyElements);
    matrix_[i]->apply(X, *curY); // This gives a result that is off-diagonal-matrix*X
  } //end for

  RCP<VEC> tempVec = rcp(new VEC(invDiagonal_->getMap()));
  tempVec->reciprocal(*invDiagonal_);

  Y.elementWiseMultiply(1.0,*tempVec, X, 1.0); // Add diagonal contribution

} //end Apply
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
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

  if (verbose)
  {
    std::cout << "A11 self-check residual = " << resid << std::endl;
  } //end if

  TEUCHOS_TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n");

} //end Check
#if LINSOLVE_PREC == 0
// Use float
template class dft_PolyA11_Tpetra_Operator<float, int, int>;
#elif LINSOLVE_PREC == 1
// Use double
template class dft_PolyA11_Tpetra_Operator<double, int, int>;
#if MIXED_PREC == 1
template class dft_PolyA11_Tpetra_Operator<float, int, int>;
#endif
#elif LINSOLVE_PREC == 2
// Use quad double
template class dft_PolyA11_Tpetra_Operator<qd_real, int, int>;
#if MIXED_PREC == 1
template class dft_PolyA11_Tpetra_Operator<dd_real, int, int>;
#endif
#elif LINSOLVE_PREC == 3
// Use double double
template class dft_PolyA11_Tpetra_Operator<dd_real, int, int>;
#if MIXED_PREC == 1
template class dft_PolyA11_Tpetra_Operator<double, int, int>;
#endif
#endif
