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
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA11_Coulomb_Tpetra_Operator
(RCP<const MAP > & ownedMap, RCP<const MAP > & block1Map, 
 RCP<const MAP > & allGMap, RCP<const MAP > & poissonMap, 
 RCP<ParameterList> parameterList) 
  : dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>(ownedMap, allGMap),
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
  poissonMatrix_->setObjectLabel("PolyA11Coulomb::poissonMatrix");
  return;
} //end constructor
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_PolyA11_Coulomb_Tpetra_Operator
() 
{
  return;
} //end destructor
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
() 
{
  
  TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n"); 
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) 
  {
    for (LocalOrdinal i=0; i<numBlocks_; i++)
    {
      matrix_[i]->setAllToScalar(0.0);
    } //end for
    poissonMatrix_->setAllToScalar(0.0);
  } //end if 
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value) 
{
  Array<GlobalOrdinal> cols(1);
  cols[0] = colGID;
  Array<Scalar> vals(1);
  vals[0] = value;
  if (ownedPhysicsID >= numBlocks_) //insert it into Poisson part
  {
    if (firstTime_) 
    {
      if (rowGID != curPoissonRow_) 
      {
	insertPoissonRow();
	curPoissonRow_ = rowGID;
	curPoissonOwnedNode_ = ownedNode;
      } //end if
      curPoissonRowValues_[colGID] += value;
    } //end if
    else 
    {
      poissonMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
    } //end else
  } //end if
  else //insert it into G part
  {
    if (rowGID!=colGID) 
    {
      value = -value; // negate off-diagonal values to simplify kernel calls
    } //end if

    if (firstTime_) 
    {
      if (rowGID!=curRow_) 
      { 
      P11TO::insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
      curOwnedPhysicsID_ = ownedPhysicsID;
      curOwnedNode_ = ownedNode;
      } //end if
    curRowValues_[colGID] += value;
    } //end if
    else
    {
      matrix_[ownedPhysicsID]->sumIntoGlobalValues(ownedNode, cols, vals);
    } //end else
  } //end else
} //end insertMatrixValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertPoissonRow
() 
{
  if (curPoissonRowValues_.empty()) return;
  size_t numEntries = curPoissonRowValues_.size();
  if (numEntries>indices_.size()) {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  }
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos=curPoissonRowValues_.begin(); pos!=curPoissonRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  poissonMatrix_->insertGlobalValues(curPoissonRow_, indices_, values_);
  curPoissonRowValues_.clear();
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
() 
{
  if (isLinearProblemSet_) 
  {
    return; // nothing to do
  } //end if
  
  if (firstTime_) 
  { 
    P11TO::insertRow();
    // Dump any remaining entries
    insertPoissonRow();
  } //end if
  for (LocalOrdinal i=0; i<numBlocks_; i++) 
  {
    matrix_[i]->fillComplete(allGMap_, ownedMap_);
    //TEST_FOR_EXCEPT(!matrix_[i]->LowerTriangular());
  } //end for
  poissonMatrix_->fillComplete(poissonMap_, poissonMap_); 

  for (LocalOrdinal i = 0; i < poissonMap_->getNodeNumElements(); i++) 
  {
    GlobalOrdinal row = poissonMatrix_->getRowMap()->getGlobalElement(i);
    Array<Scalar> values(1);
    values[0] =10.0e-12;
    Array<GlobalOrdinal> indices(1);
    indices[0] = row;
    poissonMatrix_->sumIntoGlobalValues(row, indices, values);
  } //end for

  isLinearProblemSet_ = true;
  firstTime_ = false;
} //end finalizeProblemValues
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const 
{
  TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap())); 
  TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());

#ifdef KDEBUG
  printf("\n\n\n\ndft_PolyA11_Coulomb_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  size_t NumVectors = Y.getNumVectors();
  size_t numMyElements = ownedMap_->getNodeNumElements();

  Y=X; // We can safely do this

  RCP<MV> Y2 = Y.offsetViewNonConst(poissonMap_, numMyElements*numBlocks_);
  RCP<MV> Y1 = Y.offsetViewNonConst(allGMap_, 0);
  RCP<MV> Y1tmp = Y.offsetViewNonConst(ownedMap_, 0);
  // Start Y1tmp to view first numNodes elements of Y1
  
  LocalOrdinal offsetAmount = 0;
  for (LocalOrdinal i=0; i< numBlocks_; i++) 
  {
    matrix_[i]->apply(*Y1, *Y1tmp);
    offsetAmount += numMyElements;
    Y1tmp = Y.offsetViewNonConst(ownedMap_, offsetAmount);
    // Reset view to next block
  } //end for

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
  RCP<LinPROB> problem = rcp(new LinPROB(poissonMatrix_, Y2, Y2));
  TEST_FOR_EXCEPT(problem->setProblem() == false);
  RCP<SolMGR> solver = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem, parameterList_));
  ReturnType ret = solver->solve();
#endif
} //end applyInverse
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const 
{
  TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
  size_t NumVectors = Y.getNumVectors();
  size_t numMyElements = ownedMap_->getNodeNumElements();

  RCP<const MV> X1 = X.offsetView(allGMap_, 0);
  RCP<const MV> X2 = X.offsetView(poissonMap_, 0);
  RCP<MV> Y1tmp = Y.offsetViewNonConst(ownedMap_, 0);
  // Start Y1tmp to view first numNodes elements of Y1

  RCP<const MV> Xtmp = X.offsetView(ownedMap_, 0);
   // Start Xtmp to view first numNodes elements of X

  LocalOrdinal offsetValue = 0;
  for (LocalOrdinal i=0; i< numBlocks_; i++) 
  {
    matrix_[i]->apply(*X1, *Y1tmp); // This gives a result that is X - off-diagonal-matrix*X
    Y1tmp->update(-2.0, *Xtmp, 1.0); // This gives a result of -X - off-diagonal-matrix*X
    Y1tmp->scale(-1.0); // Finally negate to get the desired result
    offsetValue += numMyElements;
    Y1tmp = Y.offsetViewNonConst(ownedMap_, offsetValue);
    // Reset view to next block
    Xtmp = X.offsetView(ownedMap_, offsetValue);
    // Reset view to next block
  } //end for

  //now to apply the poissonMatrix_ to the last chunk of X
  RCP<MV > Y2 = Y.offsetViewNonConst(poissonMap_, numMyElements*numBlocks_);

  poissonMatrix_->apply(*X2, *Y2);

} //end Apply
#if LINSOLVE_PREC == 0
// Use float
template class dft_PolyA11_Coulomb_Tpetra_Operator<float, int, int>;
#elif LINSOLVE_PREC == 1
// Use double
template class dft_PolyA11_Coulomb_Tpetra_Operator<double, int, int>;
#elif LINSOLVE_PREC == 2
// Use quad double
template class dft_PolyA11_Coulomb_Tpetra_Operator<qd_real, int, int>;
#endif

