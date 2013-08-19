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
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
dft_PolyA11_Tpetra_Operator
(const RCP<const MAP > & ownedMap, const RCP<const MAP > & block1Map,
 RCP<ParameterList> parameterList)
  : ownedMap_(ownedMap),
    block1Map_(block1Map),
    parameterList_(parameterList),
    numBlocks_(block1Map->getNodeNumElements()/ownedMap->getNodeNumElements()),
    nnz_(0),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1),
    curOwnedPhysicsID_(-1),
    curOwnedNode_(-1)
{
  Label_ = "dft_PolyA11_Tpetra_Operator";

  diagonal_ = rcp(new VEC(block1Map));
  invDiagonal_ = rcp(new VEC(block1Map));

  matrix_.resize(numBlocks_-1);
  matrixOperator_.resize(numBlocks_-1);
  for (LocalOrdinal i=OTLO::zero(); i<numBlocks_-1; i++)
  {
    matrix_[i] = rcp(new MAT(ownedMap, 0));
    matrixOperator_[i] = rcp(new MMOP(matrix_[i]));
    matrix_[i]->setObjectLabel("PolyA11::matrix[i]");
  }

  return;
} //end constructor
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
~dft_PolyA11_Tpetra_Operator
()
{
  return;
} //end destructor
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
void
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
initializeProblemValues
()
{
  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
  {
    for (LocalOrdinal i=OTLO::zero(); i<numBlocks_-1; i++)
    {
      matrix_[i]->resumeFill();
      matrix_[i]->setAllToScalar(STMS::zero());
    }

    diagonal_->putScalar(STS::zero());
    invDiagonal_->putScalar(STS::zero());
  }

} //end initializeProblemValues
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
void
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, GlobalOrdinal rowGID, GlobalOrdinal colGID, MatScalar value)
{

  if (rowGID==colGID)
  {
    LocalOrdinal locDiag = block1Map_->getLocalElement(colGID);
    diagonal_->sumIntoLocalValue(locDiag, value);
    return;
  }
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
      this->insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
      curOwnedPhysicsID_ = ownedPhysicsID;
      curOwnedNode_ = ownedNode;
    }
    curRowValues_[colGID] += value;
  }
  else
    matrix_[ownedPhysicsID-1]->sumIntoGlobalValues(ownedNode, Array<GlobalOrdinal>(1,colGID), Array<MatScalar>(1,value));

} //end insertMatrixValue
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
void
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
insertRow
()
{
  if (curRowValues_.empty())
  {
    return;
  }
  size_t numEntries = curRowValues_.size();
  if (numEntries>indices_.size())
  {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  }
  LocalOrdinal i=OTLO::zero();

  for (ITER pos = curRowValues_.begin(), e = curRowValues_.end(); pos != e; ++pos)
  {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }

  matrix_[curOwnedPhysicsID_-1]->insertGlobalValues(curOwnedNode_, indices_, values_);

  indices_.clear();
  values_.clear();
  curRowValues_.clear();
} //end insertRow
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
void
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
  {
    return;
  }

  if (firstTime_)
  {
    this->insertRow(); // Dump any remaining entries
  }
  RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));

  // Fill complete, and compute the total number of entries in the A11 block  
  for (LocalOrdinal i=OTLO::zero(); i<numBlocks_-1; i++)
  {
    matrix_[i]->fillComplete(block1Map_, ownedMap_, pl);
    if (firstTime_) {
      nnz_ += matrix_[i]->getGlobalNumEntries();
    }
  }
  invDiagonal_->reciprocal(*diagonal_); // Invert diagonal values for faster applyInverse() method

  isLinearProblemSet_ = true;
  firstTime_ = false;
} //end finalizeProblemValues
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
void
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
applyInverse
(const MV& X, MV& Y) const
{

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  printf("\n\n\n\ndft_PolyA11_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  Scalar ONE = STS::one();
  Scalar ZERO = STS::zero();

  size_t NumVectors = Y.getNumVectors();
  size_t numMyElements = ownedMap_->getNodeNumElements();
  RCP<MV > Ytmp = rcp(new MV(ownedMap_,NumVectors));

  Y=X; // We can safely do this

  RCP<MV > curY = Y.offsetViewNonConst(ownedMap_, 0);

  RCP<VEC> diagVec = invDiagonal_->offsetViewNonConst(ownedMap_, 0)->getVectorNonConst(0);

  curY->elementWiseMultiply(ONE, *diagVec, *curY, ZERO); // Scale Y by the first block diagonal

  // Loop over block 1 through numBlocks (indexing 0 to numBlocks-1)
  for (LocalOrdinal i=OTLO::zero(); i< numBlocks_-1; i++)
  {
    // Update views of Y and diagonal blocks
    curY = Y.offsetViewNonConst(ownedMap_, (i+1)*numMyElements);

    diagVec = invDiagonal_->offsetViewNonConst(ownedMap_, (i+1)*numMyElements)->getVectorNonConst(0);

    matrixOperator_[i]->apply(Y, *Ytmp); // Multiply block lower triangular block
    curY->update(-ONE, *Ytmp, ONE); // curY = curX - Ytmp (Note that curX is in curY from initial copy Y = X)
    curY->elementWiseMultiply(ONE, *diagVec, *curY, ZERO); // Scale Y by the first block diagonal
  }
} //end applyInverse
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
void
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
#endif

  size_t numMyElements = ownedMap_->getNodeNumElements();

  RCP<MV > curY = Y.offsetViewNonConst(ownedMap_, 0);

  for (LocalOrdinal i=OTLO::zero(); i< numBlocks_-1; i++) {
    curY = Y.offsetViewNonConst(ownedMap_, (i+1)*numMyElements);
    matrixOperator_[i]->apply(X, *curY); // This gives a result that is off-diagonal-matrix*X
  }

  Y.elementWiseMultiply(STS::one(),*diagonal_, X, STS::one()); // Add diagonal contribution

} //end Apply
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalSparseOps>
void
dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node,LocalSparseOps>::
Check
(bool verbose) const
{
  RCP<VEC > x = rcp(new VEC(getDomainMap()));
  RCP<VEC > b = rcp(new VEC(getRangeMap()));

  x->randomize(); // Fill x with random numbers
  apply(*x, *b); // Forward operation
  applyInverse(*b, *b); // Reverse operation
  b->update(-STS::one(), *x, STS::one()); // Should be zero

  Scalar absResid = b->norm2();
  Scalar normX = x->norm2();
  Scalar resid = absResid / normX;

  if (verbose)
  {
    std::cout << "A11 self-check residual = " << resid << std::endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n");

} //end Check

TRAMONTO_INST_HELPER(dft_PolyA11_Tpetra_Operator)
