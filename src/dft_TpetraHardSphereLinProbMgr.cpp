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
template <class Scalar, class MatrixType>
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
dft_HardSphereLinProbMgr
(LocalOrdinal numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm, RCP<Node> node, bool formSchurMatrix, bool debug)
  : dft_BasicLinProbMgr<Scalar,MatrixType>
    (numUnknownsPerNode, parameterList, comm, node),
    isA22Diagonal_(false),
    formSchurMatrix_(formSchurMatrix),
    debug_(debug),
    curRowA12_(-1),
    curRowA21_(-1)
{
#if VERB_LEVEL > 0
  printf("\n\n\nCreated a HardSphereLinProbMgr.\n\n\n");
#endif
  return;
}

//=============================================================================
template <class Scalar, class MatrixType>
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
~dft_HardSphereLinProbMgr
()
{
   return;
}

//=============================================================================
template <class Scalar, class MatrixType>
void
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
finalizeBlockStructure
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(isBlockStructureSet_, std::runtime_error, "Block structure already set.\n");
#endif
  // Create importer to map from owned nodes to box nodes

  ownedToBoxImporter_ = Teuchos::rcp(new IMP(ownedMap_, boxMap_));

#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION((numGlobalNodes_==0 ||
			      numGlobalBoxNodes_==0 ||
			      indNonLocalEquations_.size()==0 ||
			      depNonLocalEquations_.size()<0  ||
			      densityEquations_.size()==0), std::logic_error,
			     "One or more set methods not called.");
#endif

  //Not checking if poissonEquations_.Length()==0 because don't HAVE to have Poisson equations
  //Not checking if gInvEquations_.Length()==0 because we don't have to have G inv equations

  // Fill physics ordering vector with the concatenated contents of the IDs for all physics types
  // Load Schur block mappings
  physicsOrdering_.clear();
  physicsIdToSchurBlockId_.resize(numUnknownsPerNode_);

  for (LocalOrdinal i=OTLO::zero(); i<indNonLocalEquations_.size(); ++i) {
    physicsOrdering_.append(indNonLocalEquations_[i]);
    physicsIdToSchurBlockId_[indNonLocalEquations_[i]] = 1;
  }
  for (LocalOrdinal i=OTLO::zero(); i<depNonLocalEquations_.size(); ++i) {
    physicsOrdering_.append(depNonLocalEquations_[i]);
    physicsIdToSchurBlockId_[depNonLocalEquations_[i]] = 1;
  }
  for (LocalOrdinal i=OTLO::zero(); i<densityEquations_.size(); ++i) {
    physicsOrdering_.append(densityEquations_[i]);
    physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
  }

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=OTLO::zero(); i<physicsOrdering_.size(); ++i) {
    solverOrdering_[physicsOrdering_[i]]=i;
  }

  // Special setup for coarsened nodes
  // Build int vector of with 1 if node is coarsened
  ownedNodeIsCoarsened_ = rcp(new VEC(ownedMap_)); // Assume no nodes are coarsened
  boxNodeIsCoarsened_ = rcp(new VEC(boxMap_)); // Assume no nodes are coarsened

  if (numGlobalCoarsenedNodes_>0) {
    setA22BlockIsDiagonal(false); // A22 block is not diagonal when some nodes are coarsened
    for (LocalOrdinal i=OTLO::zero(); i<numCoarsenedNodes_; ++i)
      ownedNodeIsCoarsened_->replaceLocalValue(ownedMap_->getLocalElement(coarsenedNodesMap_->getGlobalElement(i)), STS::one());
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
  LocalOrdinal k=OTLO::zero();
  LocalOrdinal k1 = (numOwnedNodes_ - numCoarsenedNodes_) * numUnknownsPerNode_; // starting point for coarsened variables
  for (LocalOrdinal i=OTLO::zero(); i<numUnknownsPerNode_; ++i) {
    LocalOrdinal ii=physicsOrdering_[i];
    if (numCoarsenedNodes_==0) {
      for (LocalOrdinal j=OTLO::zero(); j<numOwnedNodes_; ++j)
	globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
    }
    else {
      for (LocalOrdinal j=OTLO::zero(); j<numOwnedNodes_; ++j) {
	LocalOrdinal curGID = GIDs[j];
	if (coarsenedNodesMap_->isNodeGlobalElement(curGID))
	  globalGIDList[k1++] = ii*numGlobalNodes_ + GIDs[j];
	else
	  globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
      }
    }
  }

  Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  globalRowMap_ = rcp(new MAP(INVALID, globalGIDList(0,numUnks), 0, comm_, node_));
  block1RowMap_ = rcp(new MAP(INVALID, globalGIDList(0,numUnks1), 0, comm_, node_));
  block2RowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1, numUnks2), 0, comm_, node_));
  indNonLocalRowMap_ = rcp(new MAP(INVALID, globalGIDList(0, numIndNonLocal), 0, comm_, node_));
  depNonLocalRowMap_ = rcp(new MAP(INVALID, globalGIDList(numIndNonLocal, numDepNonLocal), 0, comm_, node_));

  A11_ = rcp(new HS11TO(indNonLocalRowMap_, depNonLocalRowMap_, block1RowMap_, parameterList_));

  A12_ = rcp(new MAT(block1RowMap_, 0)); A12_->setObjectLabel("HardSphere::A12");
  A21_ = rcp(new MAT(block2RowMap_, 0)); A21_->setObjectLabel("HardSphere::A21");

  if (isA22Diagonal_)
    {
      A22Diagonal_ = rcp(new HS22TO(block2RowMap_));
      A22DiagonalPrecond_ = rcp(new INVOP(A22Diagonal_));
    }
  else
    {
      A22Matrix_ = rcp(new A22MTO(block2RowMap_, parameterList_));
      A22MatrixPrecond_ = rcp(new INVOP(A22Matrix_));
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

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;
}

//=============================================================================
template <class Scalar, class MatrixType>
void
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
initializeProblemValues
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(!isBlockStructureSet_, std::logic_error,
			     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
  TEUCHOS_TEST_FOR_EXCEPTION(!isGraphStructureSet_, std::logic_error,
			     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
#endif
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    A12Static_->resumeFill();
    A21Static_->resumeFill();
    A12Static_->setAllToScalar(STMS::zero());
    A21Static_->setAllToScalar(STMS::zero());
    globalRhs_->putScalar(STS::zero());
    globalLhs_->putScalar(STS::zero());
    if (debug_) {
      globalMatrix_->setAllToScalar(STMS::zero());
    }
  }

  A11_->initializeProblemValues();
  if (isA22Diagonal_)
    A22Diagonal_->initializeProblemValues();
  else
    A22Matrix_->initializeProblemValues();

}

//=============================================================================
template <class Scalar, class MatrixType>
void dft_HardSphereLinProbMgr<Scalar,MatrixType>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, MatScalar value) {

  ArrayRCP<const Scalar> ownedNodeIsCoarsenedValues = ownedNodeIsCoarsened_->get1dView();
  ArrayRCP<const Scalar> boxNodeIsCoarsenedValues = boxNodeIsCoarsened_->get1dView();
  LocalOrdinal schurBlockRow = (physicsIdToSchurBlockId_[ownedPhysicsID]==1 &&
				ownedNodeIsCoarsenedValues[ownedNode]==0) ? 1: 2;
  LocalOrdinal schurBlockCol = (physicsIdToSchurBlockId_[boxPhysicsID]==1 &&
				boxNodeIsCoarsenedValues[boxNode]==0) ? 1: 2;
  GlobalOrdinal rowGID = this->ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = this->boxToSolverGID(boxPhysicsID, boxNode);

  LocalOrdinal schurBlockNumber = 10 * schurBlockRow + schurBlockCol;

  switch (schurBlockNumber)
  {
  case 11:
    // A11 block
    A11_->insertMatrixValue(rowGID, colGID, value);
    break;
  case 22:
    // A22 block
    if (isA22Diagonal_ && numCoarsenedNodes_==0)
      A22Diagonal_->insertMatrixValue(rowGID, colGID, value);
    else
     A22Matrix_->insertMatrixValue(rowGID, colGID, value);
    break;
  case 21:
    // A21 block
    if (firstTime_) {
      if (rowGID!=curRowA21_) {
	// Insert the current row values into the matrix and move on to the next row
	this->insertRowA21();
	curRowA21_=rowGID;
      }
      curRowValuesA21_[colGID] += value;
    }
    else
      A21Static_->sumIntoGlobalValues(rowGID, Array<GlobalOrdinal>(1,colGID), Array<MatScalar>(1,value));
    break;
  case 12:
    // A12 block
    if (firstTime_) {
      if (rowGID!=curRowA12_) {
	// Insert the current row values into the matrix and move on to the next row
	this->insertRowA12();
	curRowA12_=rowGID;
      }
      curRowValuesA12_[colGID] += value;
    }
    else
      A12Static_->sumIntoGlobalValues(rowGID, Array<GlobalOrdinal>(1,colGID), Array<MatScalar>(1,value));
    break;
  }

  if (debug_) {
    BLPM::insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
  }

}

//=============================================================================
template <class Scalar, class MatrixType>
void
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
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
  LocalOrdinal i=OTLO::zero();

  for (ITER pos = curRowValuesA12_.begin(), e = curRowValuesA12_.end(); pos != e; ++pos) {
    indicesA12_[i] = pos->first;
    valuesA12_[i++] = pos->second;
  }
  A12_->insertGlobalValues(curRowA12_, indicesA12_, valuesA12_);

  indicesA12_.clear();
  valuesA12_.clear();
  curRowValuesA12_.clear();

}

//=============================================================================
template <class Scalar, class MatrixType>
void
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
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
  LocalOrdinal i=OTLO::zero();

  for (ITER pos = curRowValuesA21_.begin(), e = curRowValuesA21_.end(); pos != e; ++pos) {
    indicesA21_[i] = pos->first;
    valuesA21_[i++] = pos->second;
  }
  A21_->insertGlobalValues(curRowA21_, indicesA21_, valuesA21_);

  indicesA21_.clear();
  valuesA21_.clear();
  curRowValuesA21_.clear();

}

//=============================================================================
template <class Scalar, class MatrixType>
void
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
finalizeProblemValues
()
{

  if (isLinearProblemSet_) {
    return; // nothing to do
  }

  if (firstTime_) {
    this->insertRowA12(); // Dump any remaining entries
    this->insertRowA21(); // Dump any remaining entries

    RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
    pl->set( "Preserve Local Graph", true );
    if(!A12_->isFillComplete()){
      A12_->fillComplete(block2RowMap_,block1RowMap_,pl);
    }
    if(!A21_->isFillComplete()) {
      A21_->fillComplete(block1RowMap_,block2RowMap_,pl);
    }

    ArrayRCP<size_t> numEntriesPerRowA12(A12_->getRowMap()->getNodeNumElements());
    for (LocalOrdinal i = OTLO::zero(); i < A12_->getRowMap()->getNodeNumElements(); ++i) {
      numEntriesPerRowA12[i] = A12_->getNumEntriesInLocalRow( i );
    }
    ArrayRCP<size_t> numEntriesPerRowA21(A21_->getRowMap()->getNodeNumElements());
    for (LocalOrdinal i = OTLO::zero(); i < A21_->getRowMap()->getNodeNumElements(); ++i) {
      numEntriesPerRowA21[i] = A21_->getNumEntriesInLocalRow( i );
    }

    A12Graph_ = rcp(new GRAPH(A12_->getRowMap(), A12_->getColMap(), numEntriesPerRowA12, Tpetra::StaticProfile));
    A21Graph_ = rcp(new GRAPH(A21_->getRowMap(), A21_->getColMap(), numEntriesPerRowA21, Tpetra::StaticProfile));
    for (LocalOrdinal i = OTLO::zero(); i < A12_->getRowMap()->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      A12_->getLocalRowView( i, indices, values );
      A12Graph_->insertLocalIndices( i, indices );
    }
    for (LocalOrdinal i = OTLO::zero(); i < A21_->getRowMap()->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      A21_->getLocalRowView( i, indices, values );
      A21Graph_->insertLocalIndices( i, indices );
    }
    A12Graph_->fillComplete(block2RowMap_,block1RowMap_);
    A21Graph_->fillComplete(block1RowMap_,block2RowMap_);
    A12Static_ = rcp(new MAT(A12Graph_));
    A21Static_ = rcp(new MAT(A21Graph_));
    A12Static_->setAllToScalar(STMS::zero());
    A21Static_->setAllToScalar(STMS::zero());

    for (LocalOrdinal i = OTLO::zero(); i < A12_->getRowMap()->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      A12_->getLocalRowView( i, indices, values );
      A12Static_->sumIntoLocalValues( i, indices(), values() );
    }
    A12Static_->fillComplete(block2RowMap_,block1RowMap_,pl);
    for (LocalOrdinal i = OTLO::zero(); i < A21_->getRowMap()->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      A21_->getLocalRowView( i, indices, values );
      A21Static_->sumIntoLocalValues( i, indices(), values() );
    }
    A21Static_->fillComplete(block1RowMap_,block2RowMap_,pl);

    if (isA22Diagonal_)
      schurOperator_ = rcp(new ScTO(A11_, A12Static_, A21Static_, A22Diagonal_));
    else
      schurOperator_ = rcp(new ScTO(A11_, A12Static_, A21Static_, A22Matrix_));

    if (debug_) {
      globalMatrix_->fillComplete();
    }

  }
  RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
  if(!A12Static_->isFillComplete()){
    A12Static_->fillComplete(block2RowMap_,block1RowMap_,pl);
  }
  if(!A21Static_->isFillComplete()) {
    A21Static_->fillComplete(block1RowMap_,block2RowMap_,pl);
  }
  //std::cout << *A12_ << endl
  //          << *A21_ << endl;

  A11_->finalizeProblemValues();
  if (isA22Diagonal_)
    A22Diagonal_->finalizeProblemValues();
  else
    A22Matrix_->finalizeProblemValues();

  // Compute the dimension and total number of entries of the matrix
  if (firstTime_) {
    GlobalOrdinal nnz = 0;
    if (isA22Diagonal_)
    {
      nnz =  A22Diagonal_->getNumEntries();
    }
    else
    {
      nnz =  A22Matrix_->getNumEntries();
    }
    GlobalOrdinal dim = globalRowMap_->getGlobalNumElements();
    nnz += A11_->getNumEntries();
    nnz += A12Static_->getGlobalNumEntries();
    nnz += A21Static_->getGlobalNumEntries();

#if VERB_LEVEL > 0
    printf("\n\nGlobal matrix has %d rows and %d nonzeros..\n\n", dim, nnz);
#endif
  }

  //  Check(true);
  isLinearProblemSet_ = true;
  firstTime_ = false;

}

//=============================================================================
template <class Scalar, class MatrixType>
void
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
setupSolver
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error,
			     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");
#endif

  schurOperator_->ComputeRHS(*rhs1_, *rhs2_, *rhsSchur_);

  problem_ = rcp(new LinPROB(schurOperator_, lhs2_, rhsSchur_));

  if (formSchurMatrix_) {// We have S explicitly available, so let's use it
    if (isA22Diagonal_)
      schurOperator_->SetSchurComponents(A11_->getA11invMatrix(), A22Diagonal_->getA22Matrix());
    else
      schurOperator_->SetSchurComponents(A11_->getA11invMatrix(), A22Matrix_->getA22Matrix());
    schurComplementMatrixOperator_ = rcp(new MMOP(schurOperator_->getSchurComplement()));
    problem_->setOperator(schurComplementMatrixOperator_);
  }
  int precond  = parameterList_->template get<int>( "Precond" );
  if (precond != AZ_none) {
    if (isA22Diagonal_)
      problem_->setLeftPrec(A22DiagonalPrecond_);
    else
      problem_->setLeftPrec(A22MatrixPrecond_);
  }

  TEUCHOS_TEST_FOR_EXCEPT(problem_->setProblem() == false);
  RCP<ParameterList> belosList = rcp(new ParameterList(parameterList_->sublist("belosList")));
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem_, belosList));

}

//=============================================================================
template <class Scalar, class MatrixType>
void dft_HardSphereLinProbMgr<Scalar,MatrixType>::
solve
()
{

  // Solve the Schur complement system
  try {
    solver_->solve();
  }
  catch (Belos::StatusTestError& e) {
    std::cout << "Belos failed to solve the linear problem! Belos threw exception "
	      << e.what() << std::endl;
  }

  // Compute the rest of the solution
  schurOperator_->ComputeX1(*rhs1_, *lhs2_, *lhs1_);

  if (debug_) {
    RCP<VEC> tmpRhs = rcp(new VEC(globalRowMap_));
    RCP<VEC> tmprhs1 = tmpRhs->offsetViewNonConst(block1RowMap_, 0)->getVectorNonConst(0);
    RCP<VEC> tmprhs2 = tmpRhs->offsetViewNonConst(block2RowMap_, block1RowMap_->getNodeNumElements())->getVectorNonConst(0);

    schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *tmprhs1, *tmprhs2);

    tmpRhs->update(-STS::one(), *globalRhs_, STS::one());
    Scalar resid = STS::zero();
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

}

//=============================================================================
template <class Scalar, class MatrixType>
ArrayRCP<ArrayRCP<Scalar> >
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
applyMatrix
(const ArrayView<const ArrayView<const Scalar> >& x) const
{
  this->setLhs(x);

  schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *rhs1_, *rhs2_);

  return (getRhs());
}

//=============================================================================
template <class Scalar, class MatrixType>
void dft_HardSphereLinProbMgr<Scalar,MatrixType>::
Check
(bool verbose) const
{
  A11_->Check(verbose);
  if (isA22Diagonal_) {
    A22Diagonal_->Check(verbose);
  }

}

//=============================================================================
template <class Scalar, class MatrixType>
void
dft_HardSphereLinProbMgr<Scalar,MatrixType>::
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

TRAMONTO_INST_HELPER(dft_HardSphereLinProbMgr)
