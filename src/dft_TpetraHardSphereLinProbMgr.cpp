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

#include "dft_TpetraHardSphereLinProbMgr.hpp"

//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_HardSphereLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_HardSphereLinProbMgr
(LocalOrdinal numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm, bool formSchurMatrix, bool debug)
  : dft_BasicLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>
    (numUnknownsPerNode, parameterList, comm),
    isA22Diagonal_(false),
    formSchurMatrix_(formSchurMatrix),
    debug_(debug),
    curRowA12_(-1),
    curRowA21_(-1) {
  return;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
~dft_HardSphereLinProbMgr
()
{
   return;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
finalizeBlockStructure
()
{
  TEUCHOS_TEST_FOR_EXCEPTION(isBlockStructureSet_, std::runtime_error, "Block structure already set.\n");

  // Create importer to map from owned nodes to box nodes

  ownedToBoxImporter_ = Teuchos::rcp(new IMP(ownedMap_, boxMap_));

  TEUCHOS_TEST_FOR_EXCEPTION((numGlobalNodes_==0 ||
			      numGlobalBoxNodes_==0 ||
			      indNonLocalEquations_.size()==0 ||
			      depNonLocalEquations_.size()<0  ||
			      densityEquations_.size()==0), std::logic_error,
			     "One or more set methods not called.");
  //Not checking if poissonEquations_.Length()==0 because don't HAVE to have Poisson equations
  //Not checking if gInvEquations_.Length()==0 because we don't have to have G inv equations

  // Fill physics ordering vector with the concatenated contents of the IDs for all physics types
  // Load Schur block mappings
  physicsOrdering_.clear();
  physicsIdToSchurBlockId_.resize(numUnknownsPerNode_);

  for (LocalOrdinal i=0; i<indNonLocalEquations_.size(); i++) {
    physicsOrdering_.append(indNonLocalEquations_[i]);
    physicsIdToSchurBlockId_[indNonLocalEquations_[i]] = 1;
  }
  for (LocalOrdinal i=0; i<depNonLocalEquations_.size(); i++) {
    physicsOrdering_.append(depNonLocalEquations_[i]);
    physicsIdToSchurBlockId_[depNonLocalEquations_[i]] = 1;
  }
  for (LocalOrdinal i=0; i<densityEquations_.size(); i++) {
    physicsOrdering_.append(densityEquations_[i]);
    physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
  }

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++) {
    solverOrdering_[physicsOrdering_[i]]=i;
  }

  // Special setup for coarsened nodes
  // Build int vector of with 1 if node is coarsened
  ownedNodeIsCoarsened_ = rcp(new VEC(ownedMap_)); // Assume no nodes are coarsened
  boxNodeIsCoarsened_ = rcp(new VEC(boxMap_)); // Assume no nodes are coarsened

  if (numGlobalCoarsenedNodes_>0) {
    setA22BlockIsDiagonal(false); // A22 block is not diagonal when some nodes are coarsened
    for (LocalOrdinal i=0; i<numCoarsenedNodes_; i++)
      ownedNodeIsCoarsened_->replaceLocalValue(ownedMap_->getLocalElement(coarsenedNodesMap_->getGlobalElement(i)), 1.0);
    boxNodeIsCoarsened_->doImport(*ownedNodeIsCoarsened_, *ownedToBoxImporter_, Tpetra::INSERT); // Now each processor knows which of its box nodes is coarsened
  }
  const size_t numUnks = numOwnedNodes_*numUnknownsPerNode_;
  const size_t numRealNodes = numOwnedNodes_ - numCoarsenedNodes_;
  const size_t numUnks1 = numRealNodes*(indNonLocalEquations_.size()+depNonLocalEquations_.size());
  const size_t numUnks2 = numRealNodes*(densityEquations_.size()) + numCoarsenedNodes_ * numUnknownsPerNode_;
  assert(numUnks==(numUnks1+numUnks2));  // Sanity test
  const size_t numIndNonLocal = numRealNodes*(indNonLocalEquations_.size());
  const size_t numDepNonLocal = numRealNodes*(depNonLocalEquations_.size());
  const int numDensity = numRealNodes*(densityEquations_.size());
  Array<GlobalOrdinal> globalGIDList(numUnks);

  ArrayView<const GlobalOrdinal> GIDs = ownedMap_->getNodeElementList();
  LocalOrdinal k=0;
  LocalOrdinal k1 = (numOwnedNodes_ - numCoarsenedNodes_) * numUnknownsPerNode_; // starting point for coarsened variables
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    LocalOrdinal ii=physicsOrdering_[i];
    if (numCoarsenedNodes_==0) {
      for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
	globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
    }
    else {
      for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
	LocalOrdinal curGID = GIDs[j];
	if (coarsenedNodesMap_->isNodeGlobalElement(curGID))
	  globalGIDList[k1++] = ii*numGlobalNodes_ + GIDs[j];
	else
	  globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
      }
    }
  }

  globalRowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList(0,numUnks), 0, comm_));
  block1RowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList(0,numUnks1), 0, comm_));
  block2RowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList(numUnks1, numUnks2), 0, comm_));
  indNonLocalRowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList(0, numIndNonLocal), 0, comm_));
  depNonLocalRowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList(numIndNonLocal, numDepNonLocal), 0, comm_));

  /*
    std::cout << " Global Row Map" << *globalRowMap_ << std::endl
    << " Block 1     Row Map " << *block1RowMap_ << std::endl
    << " Block 2     Row Map " << *block2RowMap_ << std::endl
    << " DepNonLocal Row Map " << *depNonLocalRowMap_ << std::endl;
  */

  A11_ = rcp(new HS11TO(indNonLocalRowMap_, depNonLocalRowMap_, block1RowMap_));
  A12_ = rcp(new MAT(block1RowMap_, 0)); A12_->setObjectLabel("HardSphere::A12");
  A21_ = rcp(new MAT(block2RowMap_, 0)); A21_->setObjectLabel("HardSphere::A21");
  if (isA22Diagonal_) {
    A22Diagonal_ = rcp(new HS22TO(block2RowMap_));
  }
  else {
    A22Matrix_ = rcp(new A22MTO(block2RowMap_, parameterList_));
  }
  if (debug_) {
    globalMatrix_ = rcp(new MAT(globalRowMap_, 0));
    globalMatrix_->setObjectLabel("HardSphere::globalMatrix");
    }
  else
    globalMatrix_ = Teuchos::null; // not used by this solver

  globalRhs_ = rcp(new VEC(globalRowMap_));
  globalLhs_ = rcp(new VEC(globalRowMap_));

  rhs1_ = globalRhs_->offsetViewNonConst(block1RowMap_, 0)->getVectorNonConst(0);
  rhs2_ = globalRhs_->offsetViewNonConst(block2RowMap_, numUnks1)->getVectorNonConst(0);
  rhsSchur_ = rcp(new VEC(*rhs2_));
  lhs1_ = globalLhs_->offsetViewNonConst(block1RowMap_, 0)->getVectorNonConst(0);
  lhs2_ = globalLhs_->offsetViewNonConst(block2RowMap_, numUnks1)->getVectorNonConst(0);

  if (isA22Diagonal_)
    schurOperator_ = rcp(new ScTO(A11_, A12_, A21_, A22Diagonal_));
  else
    schurOperator_ = rcp(new ScTO(A11_, A12_, A21_, A22Matrix_));

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
initializeProblemValues
()
{

  TEUCHOS_TEST_FOR_EXCEPTION(!isBlockStructureSet_, std::logic_error,
			     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
  TEUCHOS_TEST_FOR_EXCEPTION(!isGraphStructureSet_, std::logic_error,
			     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    A12_->resumeFill();
    A21_->resumeFill();
    A12_->setAllToScalar(0.0);
    A21_->setAllToScalar(0.0);
    globalRhs_->putScalar(0.0);
    globalLhs_->putScalar(0.0);
    if (debug_) {
      globalMatrix_->setAllToScalar(0.0);
    }
  }

  A11_->initializeProblemValues();
  if (isA22Diagonal_)
    A22Diagonal_->initializeProblemValues();
  else
    A22Matrix_->initializeProblemValues();

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, Scalar value) {

  ArrayRCP<const Scalar> ownedNodeIsCoarsenedValues = ownedNodeIsCoarsened_->get1dView();
  ArrayRCP<const Scalar> boxNodeIsCoarsenedValues = boxNodeIsCoarsened_->get1dView();
  bool schurBlockRow1 = (physicsIdToSchurBlockId_[ownedPhysicsID]==1 &&
			 ownedNodeIsCoarsenedValues[ownedNode]==0);
  bool schurBlockCol1 = (physicsIdToSchurBlockId_[boxPhysicsID]==1 &&
			 boxNodeIsCoarsenedValues[boxNode]==0);
  GlobalOrdinal rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = boxToSolverGID(boxPhysicsID, boxNode);

  Array<GlobalOrdinal> cols(1);
  cols[0] = colGID;
  Array<Scalar> vals(1);
  vals[0] = value;

  if (schurBlockRow1 && schurBlockCol1) { // A11 block
    A11_->insertMatrixValue(rowGID, colGID, value);
  }
  else if (!schurBlockRow1 && !schurBlockCol1) { // A22 block
    if (isA22Diagonal_ && numCoarsenedNodes_==0)
      A22Diagonal_->insertMatrixValue(rowGID, colGID, value);
    else
     A22Matrix_->insertMatrixValue(rowGID, colGID, value);
  }
  else if (!schurBlockRow1 && schurBlockCol1) { // A21 block
    if (firstTime_) {
      if (rowGID!=curRowA21_) {
	insertRowA21();  // Dump the current contents of curRowValues_ into matrix and clear map
	curRowA21_=rowGID;
      }
      curRowValuesA21_[colGID] += value;
    }
    else
      A21_->sumIntoGlobalValues(rowGID, cols, vals);
  }
  else { // A12 block
    if (firstTime_) {
      if (rowGID!=curRowA12_) {
	insertRowA12();  // Dump the current contents of curRowValues_ into matrix and clear map
	curRowA12_=rowGID;
      }
      curRowValuesA12_[colGID] += value;
    }
    else
      A12_->sumIntoGlobalValues(rowGID, cols, vals);
  }

  if (debug_) {
    BLPM::insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
  }

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
insertRowA12
()
{

  if (curRowValuesA12_.empty()) {
    return;
  }
  size_t numEntries = curRowValuesA12_.size();
  if (numEntries>indicesA12_.size()) {
    indicesA12_.resize(numEntries);
    valuesA12_.resize(numEntries);
  }
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos = curRowValuesA12_.begin(); pos != curRowValuesA12_.end(); ++pos) {
    indicesA12_[i] = pos->first;
    valuesA12_[i++] = pos->second;
  }
  A12_->insertGlobalValues(curRowA12_, indicesA12_, valuesA12_);

  indicesA12_.clear();
  valuesA12_.clear();
  curRowValuesA12_.clear();

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
insertRowA21
()
{

  if (curRowValuesA21_.empty()) {
    return;
  }
  size_t numEntries = curRowValuesA21_.size();
  if (numEntries>indicesA21_.size()) {
    indicesA21_.resize(numEntries);
    valuesA21_.resize(numEntries);
  }
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos = curRowValuesA21_.begin(); pos != curRowValuesA21_.end(); ++pos) {
    indicesA21_[i] = pos->first;
    valuesA21_[i++] = pos->second;
  }
  A21_->insertGlobalValues(curRowA21_, indicesA21_, valuesA21_);

  indicesA21_.clear();
  valuesA21_.clear();
  curRowValuesA21_.clear();

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
finalizeProblemValues
()
{

  if (isLinearProblemSet_) {
    return; // nothing to do
  }

  if (firstTime_) {
    insertRowA12(); // Dump any remaining entries
    insertRowA21(); // Dump any remaining entries

    if (debug_) {
      globalMatrix_->fillComplete();
    }
  }
  //std::cout << *A12_ << endl
  //          << *A21_ << endl;

  if(!A12_->isFillComplete()){
    A12_->fillComplete(block2RowMap_,block1RowMap_);
  }
  if(!A21_->isFillComplete()) {
    A21_->fillComplete(block1RowMap_,block2RowMap_);
  }

  A11_->finalizeProblemValues();
  if (isA22Diagonal_)
    A22Diagonal_->finalizeProblemValues();
  else
    A22Matrix_->finalizeProblemValues();

  //  Check(true);
  isLinearProblemSet_ = true;
  firstTime_ = false;

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setupSolver
()
{

  TEUCHOS_TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error,
			     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");

  schurOperator_->ComputeRHS(*rhs1_, *rhs2_, *rhsSchur_);

  ////  if (solver_ != Teuchos::null ) return(0);  //Already setup

  problem_ = rcp(new LinPROB(schurOperator_, lhs2_, rhsSchur_));

  if (formSchurMatrix_) {// We have S explicitly available, so let's use it
    if (isA22Diagonal_)
      schurOperator_->SetSchurComponents(A11_->getA11invMatrix(), A22Diagonal_->getA22Matrix());
    else
      schurOperator_->SetSchurComponents(A11_->getA11invMatrix(), A22Matrix_->getA22Matrix());
    problem_->setOperator(schurOperator_->getSchurComplement());
  }
  if (isA22Diagonal_)
    problem_->setLeftPrec(A22Diagonal_);
  else
    problem_->setLeftPrec(A22Matrix_);

  TEUCHOS_TEST_FOR_EXCEPT(problem_->setProblem() == false);
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem_, parameterList_));

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
solve
()
{

  //writeMatrix("2D.mm", "Small HardSpheremer Matrix", "Global Matrix from Small HardSpheremer Problem");
  //abort();

  //  const int * options = solver_->GetAllAztecOptions();
  //  const double * params = solver_->GetAllAztecParams();

  ReturnType ret = solver_->solve();

  schurOperator_->ComputeX1(*rhs1_, *lhs2_, *lhs1_); // Compute rest of solution

  if (debug_) {
    RCP<VEC> tmpRhs = rcp(new VEC(globalRowMap_));
    RCP<VEC> tmprhs1 = tmpRhs->offsetViewNonConst(block1RowMap_, 0)->getVectorNonConst(0);
    RCP<VEC> tmprhs2 = tmpRhs->offsetViewNonConst(block2RowMap_, block1RowMap_->getNodeNumElements())->getVectorNonConst(0);

    schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *tmprhs1, *tmprhs2);

    tmpRhs->update(-1.0, *globalRhs_, 1.0);
    Scalar resid=0.0;
    resid = tmpRhs->norm2();
    std::cout << "Global Residual for solution = " << resid << std::endl;
    bool writeMatrixNow = true;
    if (writeMatrixNow) {
      writeMatrix("A.dat", "GlobalMatrix", "GlobalMatrix");
      BLPM::writeLhs("x.dat");
      BLPM::writeRhs("b.dat");
      BLPM::writePermutation("p.dat");
    }
  }
  //std::cout << "Global RHS = " << *globalRhs_ << std::endl
  //          << "Global LHS = " << *globalLhs_ << std::endl;

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
applyMatrix
(const ArrayView<const ArrayView<const Scalar> >& x) const
{
  setLhs(x);
  schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *rhs1_, *rhs2_);
  return (getRhs());
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Check
(bool verbose) const
{
  A11_->Check(verbose);
  if (isA22Diagonal_) {
    A22Diagonal_->Check(verbose);
  }

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_HardSphereLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
writeMatrix
(const char * filename, const char * matrixName, const char * matrixDescription) const
{
  /*
  cout << "IN WRITEMATRIX ROUTINE" <<endl;
  if (debug_){
  cout << "DEBUG FLAG SEEMS TO BE TRUE"<<endl;
    return(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_, matrixName, matrixDescription));
  }
  else{
  cout << "DEBUG FLAG SEEMS TO BE false"<<endl;
    return(-1); // Not available if not in debug mode
  }
  */
}
#if LINSOLVE_PREC == 0
// Use float
template class dft_HardSphereLinProbMgr<float, int, int>;
#elif LINSOLVE_PREC == 1
// Use double
template class dft_HardSphereLinProbMgr<double, int, int>;
#elif LINSOLVE_PREC == 2
// Use quad double
template class dft_HardSphereLinProbMgr<qd_real, int, int>;
#elif LINSOLVE_PREC == 3
// Use double double
template class dft_HardSphereLinProbMgr<dd_real, int, int>;
#endif
