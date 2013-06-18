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

#include "dft_PolyA11_Coulomb_Tpetra_Operator.hpp"

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA11_Coulomb_Tpetra_Operator
(RCP<const MAP > & ownedMap, RCP<const MAP > & block1Map,
 RCP<const MAP > & allGMap, RCP<const MAP > & poissonMap,
 RCP<ParameterList> parameterList)
  : dft_PolyA11_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>(ownedMap, allGMap, parameterList),
    //dft_PolyA11_Tpetra_Operator(ownedMap, allGMap, parameterList),
    //dft_PolyA11_Tpetra_Operator(ownedMap, block1Map),
    parameterList_(parameterList),
    allGMap_(allGMap),
    poissonMap_(poissonMap),
    block1Map_(block1Map),
    curPoissonRow_(-1)
{

  Label_ = "dft_PolyA11_Coulomb_Tpetra_Operator";
  poissonMatrix_ = rcp(new MAT(poissonMap, 0));
  poissonMatrixOperator_ = rcp(new MMOP(poissonMatrix_));
  poissonMatrix_->setObjectLabel("PolyA11Coulomb::poissonMatrix");
  return;
} //end constructor
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_PolyA11_Coulomb_Tpetra_Operator
()
{
  return;
} //end destructor
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{

  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
  {
    for (LocalOrdinal i=OTLO::zero(); i<numBlocks_; i++)
    {
      matrix_[i]->setAllToScalar(STMS::zero());
    }
    poissonMatrix_->setAllToScalar(STMS::zero());
  }
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, GlobalOrdinal rowGID, GlobalOrdinal colGID, MatScalar value)
{

  if (ownedPhysicsID >= numBlocks_) //insert it into Poisson part
  {
    if (firstTime_)
    {
      if (rowGID != curPoissonRow_)
      {
	this->insertPoissonRow();
	curPoissonRow_ = rowGID;
	curPoissonOwnedNode_ = ownedNode;
      }
      curPoissonRowValues_[colGID] += value;
    }
    else
      poissonMatrix_->sumIntoGlobalValues(rowGID, Array<GlobalOrdinal>(1,colGID), Array<MatScalar>(1,value));

  }
  else //insert it into G part
  {
    if (rowGID!=colGID)
    {
      value = -value; // negate off-diagonal values to simplify kernel calls
    }

    if (firstTime_)
    {
      if (rowGID!=curRow_)
      {
      P11TO::insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
      curOwnedPhysicsID_ = ownedPhysicsID;
      curOwnedNode_ = ownedNode;
      }
      curRowValues_[colGID] += value;
    }
    else
      matrix_[ownedPhysicsID]->sumIntoGlobalValues(ownedNode, Array<GlobalOrdinal>(1,colGID), Array<MatScalar>(1,value));

  }
} //end insertMatrixValues
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertPoissonRow
()
{
  if (curPoissonRowValues_.empty()) return;
  size_t numEntries = curPoissonRowValues_.size();
  if (numEntries>indices_.size()) {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  }
  LocalOrdinal i=OTLO::zero();
  ITER pos;
  for (pos=curPoissonRowValues_.begin(); pos!=curPoissonRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  poissonMatrix_->insertGlobalValues(curPoissonRow_, indices_, values_);
  curPoissonRowValues_.clear();
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
  {
    return; // nothing to do
  }

  if (firstTime_)
  {
    P11TO::insertRow();
    // Dump any remaining entries
    this->insertPoissonRow();
  }
  RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
  for (LocalOrdinal i=OTLO::zero(); i<numBlocks_; i++)
  {
    matrix_[i]->fillComplete(allGMap_, ownedMap_, pl);
  }
  poissonMatrix_->fillComplete(poissonMap_, poissonMap_, pl);

  for (LocalOrdinal i = OTLO::zero(); i < poissonMap_->getNodeNumElements(); i++)
  {
    GlobalOrdinal row = poissonMatrix_->getRowMap()->getGlobalElement(i);
    poissonMatrix_->sumIntoGlobalValues(row, Array<GlobalOrdinal>(1,row), Array<MatScalar>(1,(MatScalar)1e-12));
  }

  problem_ = rcp(new LinPROB());
  int precond  = parameterList_->template get<int>( "Precond" );
  if (precond != AZ_none) {
    LocalOrdinal overlapLevel = 0;
    preconditioner_ = rcp(new SCHWARZ(poissonMatrix_,overlapLevel));
    preconditionerOp_ = rcp(new SCHWARZ_OP(preconditioner_));
    ParameterList ifpack2List = parameterList_->sublist("ifpack2ListA11");
    preconditioner_->setParameters(ifpack2List);
    preconditioner_->initialize();
    preconditioner_->compute();
    problem_->setLeftPrec(preconditionerOp_);
  }
  RCP<ParameterList> belosList = rcp(new ParameterList(parameterList_->sublist("belosListA11")));
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>());
  solver_->setParameters(belosList);

  isLinearProblemSet_ = true;
  firstTime_ = false;
} //end finalizeProblemValues
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const
{

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  printf("\n\n\n\ndft_PolyA11_Coulomb_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  size_t numMyElements = ownedMap_->getNodeNumElements();

  Y=X; // We can safely do this

  RCP<MV> Y2 = Y.offsetViewNonConst(poissonMap_, numMyElements*numBlocks_);
  RCP<MV> Y1 = Y.offsetViewNonConst(allGMap_, 0);
  RCP<MV> Y1tmp = Y.offsetViewNonConst(ownedMap_, 0);

  LocalOrdinal offsetAmount = 0;
  for (LocalOrdinal i=OTLO::zero(); i<numBlocks_; i++)
  {
    matrixOperator_[i]->apply(*Y1, *Y1tmp);
    offsetAmount += numMyElements;
    Y1tmp = Y.offsetViewNonConst(ownedMap_, offsetAmount);
    // Reset view to next block
  }

#ifdef SUPPORTS_STRATIMIKOS
  RCP<ThyraMV> thyraY = createMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Y2);
  RCP<ThyraOP> thyraOp = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(poissonMatrix_);
  RCP<DefaultLinearSolverBuilder> solver = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver->readParameters(out.get());

  RCP<ThyraLOWSFactory> lowsFactory = solver->createLinearSolveStrategy("");
  RCP<ThyraLOWS> lows = linearOpWithSolve<Scalar>(*lowsFactory, thyraOp);

  SolveStatus<Scalar> status = lows->solve(Thyra::NOTRANS, *thyraY, thyraY.ptr());
#else

  problem_->setOperator(poissonMatrixOperator_);
  problem_->setLHS(Y2);
  problem_->setRHS(Y2);
  TEUCHOS_TEST_FOR_EXCEPT(problem_->setProblem() == false);
  solver_->setProblem(problem_);
  solver_->solve();

#endif
} //end applyInverse
//==============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{

  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
#endif

  Scalar ONE = STS::one();

  size_t numMyElements = ownedMap_->getNodeNumElements();

  RCP<const MV> X1 = X.offsetView(allGMap_, 0);
  RCP<const MV> X2 = X.offsetView(poissonMap_, 0);
  RCP<MV> Y1tmp = Y.offsetViewNonConst(ownedMap_, 0);

  RCP<const MV> Xtmp = X.offsetView(ownedMap_, 0);

  LocalOrdinal offsetValue = 0;
  for (LocalOrdinal i=OTLO::zero(); i<numBlocks_; i++)
  {
    matrixOperator_[i]->apply(*X1, *Y1tmp); // This gives a result that is X - off-diagonal-matrix*X
    Y1tmp->update(-(ONE+ONE), *Xtmp, ONE); // This gives a result of -X - off-diagonal-matrix*X
    Y1tmp->scale(-ONE); // Finally negate to get the desired result
    offsetValue += numMyElements;
    Y1tmp = Y.offsetViewNonConst(ownedMap_, offsetValue);
    Xtmp = X.offsetView(ownedMap_, offsetValue);
    // Reset view to next block
  }

  //now to apply the poissonMatrix_ to the last chunk of X
  RCP<MV > Y2 = Y.offsetViewNonConst(poissonMap_, numMyElements*numBlocks_);

  poissonMatrixOperator_->apply(*X2, *Y2);

} //end Apply

TRAMONTO_INST_HELPER(dft_PolyA11_Coulomb_Tpetra_Operator)
