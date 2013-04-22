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

#include "dft_TpetraPolyLinProbMgr.hpp"

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
dft_PolyLinProbMgr
(LocalOrdinal numUnknownsPerNode, RCP<ParameterList> parameterList,
 RCP<const COMM> comm, RCP<Node> node, bool debug)
  : dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>
    (numUnknownsPerNode, parameterList, comm, node),
    isLinear_(false),
    debug_(debug),
    hasPoisson_(false),
    curRowA12_(-1),
    curRowA21_(-1)
{
  parameterList_->set("P_location", 0); //change
  parameterList_->set("F_location", 0); //change
#ifdef KDEBUG
  printf("\n\n\nCreated a PolyLinProbMgr.\n\n\n");
#endif
  return;
} //end constructor
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
~dft_PolyLinProbMgr
()
{
  return;
} //end destructor
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
finalizeBlockStructure
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(isBlockStructureSet_, std::runtime_error, "Block structure already set.\n");

  TEUCHOS_TEST_FOR_EXCEPTION((numGlobalNodes_==0 || numGlobalBoxNodes_==0 ||
		      gEquations_.size()==0 || cmsEquations_.size()==0 ||
		      densityEquations_.size()==0),
		      std::logic_error, "One or more set methods not called.");
#endif
  //Not checking if poissonEquations_.Length()==0 because don't HAVE to have Poisson equations
  //Not checking if gInvEquations_.Length()==0 because we don't have to have G inv equations

  // RN_20100111: Switch to the local parameter list.

  int poissonInA11 = Teuchos::getParameter<LocalOrdinal>(*parameterList_, "P_location");
  int F_location = Teuchos::getParameter<LocalOrdinal>(*parameterList_, "F_location");

  // Fill physics ordering vector with the concatenated contents of the IDs for all physics types
  // Load Schur block mappings
  physicsOrdering_.clear();
  physicsIdToSchurBlockId_.resize(numUnknownsPerNode_);
  isCmsEquation_.resize(numUnknownsPerNode_);
  isDensityEquation_.resize(numUnknownsPerNode_);
  isPoissonEquation_.resize(numUnknownsPerNode_);

  for (LocalOrdinal i=0; i<gEquations_.size(); i++)
  {
    physicsOrdering_.append(gEquations_[i]);
    physicsIdToSchurBlockId_[gEquations_[i]] = 1;
  } //end for
  for (LocalOrdinal i=0; i<gInvEquations_.size(); i++)
  {
    physicsOrdering_.append(gInvEquations_[i]);
    physicsIdToSchurBlockId_[gInvEquations_[i]] = 1;
  } //end for

  if (poissonInA11)
  {
    for (LocalOrdinal i=0; i<poissonEquations_.size(); i++)
    {
      physicsOrdering_.append(poissonEquations_[i]);
      physicsIdToSchurBlockId_[poissonEquations_[i]] = 1; //so it's in A11
      isPoissonEquation_[poissonEquations_[i]] = 1;
    } //end for
  } //end if
  else
  {
    for (LocalOrdinal i=0; i<poissonEquations_.size(); i++)
    {
      physicsOrdering_.append(poissonEquations_[i]);
      (physicsIdToSchurBlockId_)[poissonEquations_[i]] = 2; //so it's in A22
      isPoissonEquation_[poissonEquations_[i]] = 1;
    } //end for
  } //end else

  if (F_location == 1)  //F in NE
  {
    for (LocalOrdinal i=0; i<cmsEquations_.size(); i++)
    {
      physicsOrdering_.append(cmsEquations_[i]);
      physicsIdToSchurBlockId_[cmsEquations_[i]] = 2;
      isCmsEquation_[cmsEquations_[i]] = 1;
    } //end for
    for (LocalOrdinal i=0; i<densityEquations_.size(); i++)
    {
      physicsOrdering_.append(densityEquations_[i]);
      physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
      isDensityEquation_[densityEquations_[i]] = 1;
    } //end for
  } //end if
  else  //F in SW
  {
    for (LocalOrdinal i=0; i<densityEquations_.size(); i++)
    {
      physicsOrdering_.append(densityEquations_[i]);
      physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
      isDensityEquation_[densityEquations_[i]] = 1;
    } //end for
    for (LocalOrdinal i=0;  i<cmsEquations_.size(); i++)
    {
      physicsOrdering_.append(cmsEquations_[i]);
      physicsIdToSchurBlockId_[cmsEquations_[i]] = 2;
      isCmsEquation_[cmsEquations_[i]] = 1;
    } //end for
  } //end else

  // Sanity check of physics ordering
  //  BLPM::checkPhysicsOrdering();

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++)
  {
    solverOrdering_[physicsOrdering_[i]]=i;
  } //end for

  const size_t numUnks = numOwnedNodes_*numUnknownsPerNode_;
  const size_t numUnks1 = numOwnedNodes_*(gEquations_.size()+gInvEquations_.size());
  const size_t numUnks2 = numOwnedNodes_*(cmsEquations_.size()+densityEquations_.size());
  const size_t numUnksP = numOwnedNodes_*(poissonEquations_.size());
  if (numUnksP > 0)
  {
    hasPoisson_ = true;
  } //end if
  assert(numUnks==(numUnks1+numUnks2+numUnksP)); //Sanity test
  const size_t numCms = numOwnedNodes_*(cmsEquations_.size());
  const size_t numDensity = numOwnedNodes_*(densityEquations_.size());
  Array<GlobalOrdinal> globalGIDList(numUnks);

  ArrayView<const GlobalOrdinal> GIDs = ownedMap_->getNodeElementList();

  LocalOrdinal k=0;
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++)
  {
    LocalOrdinal ii=physicsOrdering_[i];
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
    } //end for
  } //end for

  Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  globalRowMap_ = rcp(new MAP(INVALID, globalGIDList(0,numUnks), 0, comm_, node_));
  if (poissonInA11)
  {
    block1RowMap_ = rcp(new MAP(INVALID, globalGIDList(0, numUnks1+numUnksP), 0, comm_, node_));
    block2RowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1+numUnksP, numUnks2), 0, comm_, node_));
  } //end if
  else
  {
    block1RowMap_ = rcp(new MAP(INVALID, globalGIDList(0, numUnks1), 0, comm_, node_));
    block2RowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1, numUnks2+numUnksP), 0, comm_, node_));
  } //end else

  if (F_location == 1) //F in NE
  {
    cmsRowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1+numUnksP, numCms), 0, comm_, node_));
    densityRowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1+numUnksP+numCms, numDensity), 0, comm_, node_));
  } //end if
  else  //F in SW
  {
    densityRowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1+numUnksP, numDensity), 0, comm_, node_));
    cmsRowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1+numUnksP+numDensity, numCms), 0, comm_, node_));
  } //end else
  if (hasPoisson_)
  {
    poissonRowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1, numUnksP), 0, comm_, node_));
    if (poissonInA11)
    {
      extraRowMap_ = rcp(new MAP(INVALID, globalGIDList(0, numUnks1), 0, comm_, node_));
      A11_ = rcp(new P11CO(ownedMap_, block1RowMap_,extraRowMap_, poissonRowMap_, parameterList_));
    } //end if
    else  //poisson in A22
    {
      extraRowMap_ = rcp(new MAP(INVALID, globalGIDList(numUnks1+numUnksP, numUnks2), 0, comm_, node_));
      A11_ = rcp(new P11TO(ownedMap_, block1RowMap_, parameterList_));
    } //end else
  } //end if
  else  //does not have Poisson equations
  {
    poissonRowMap_ = Teuchos::null;
    extraRowMap_ = Teuchos::null;
    A11_ = rcp(new P11TO(ownedMap_, block1RowMap_, parameterList_));
  } //end else

  A12_ = rcp(new MAT(block1RowMap_, 0)); A12_->setObjectLabel("PolyLinProbMgr::A12");
  A21_ = rcp(new MAT(block2RowMap_, 0)); A21_->setObjectLabel("PolyLinProbMgr::A21");

  if (hasPoisson_ && !poissonInA11)
  {
    A22_ = rcp(new P22CO(cmsRowMap_, densityRowMap_, poissonRowMap_, extraRowMap_, block2RowMap_, parameterList_));
  }
  else
  {
    A22_ = rcp(new P22TO(cmsRowMap_, densityRowMap_, block2RowMap_, parameterList_));
  }

  A22precond_ = rcp(new INVOP(A22_));

  if (debug_)
  {
    globalMatrix_ = rcp(new MAT(globalRowMap_, 0));
    globalMatrix_->setObjectLabel("PolyLinProbMgr::globalMatrix");
  } //end if
  else
  {
    globalMatrix_ = Teuchos::null; // not used by this solver
  } //end else

  globalRhs_ = rcp(new VEC(globalRowMap_));
  globalLhs_ = rcp(new VEC(globalRowMap_));

  LocalOrdinal offset2 = numUnks1;

  if (hasPoisson_ && poissonInA11)
  {
    offset2 += numUnksP;
  } //end if
  rhs1_ = globalRhs_->offsetViewNonConst(block1RowMap_, 0)->getVectorNonConst(0);
  rhs2_ = globalRhs_->offsetViewNonConst(block2RowMap_, offset2)->getVectorNonConst(0);
  rhsSchur_ = rcp(new VEC(*rhs2_));
  lhs1_ = globalLhs_->offsetViewNonConst(block1RowMap_, 0)->getVectorNonConst(0);
  lhs2_ = globalLhs_->offsetViewNonConst(block2RowMap_, offset2)->getVectorNonConst(0);

  // RN_20100113: The following is no longer needed.

  ownedToBoxImporter_ = rcp(new IMP(ownedMap_, boxMap_));

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;

} //end finalizeBlockStructure
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
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

  if (!firstTime_)
  {
    A12Static_->resumeFill();
    A21Static_->resumeFill();
    A12Static_->setAllToScalar(STMS::zero());
    A21Static_->setAllToScalar(STMS::zero());
    globalRhs_->putScalar(STS::zero());
    globalLhs_->putScalar(STS::zero());
    if (debug_)
    {
      globalMatrix_->setAllToScalar(STMS::zero());
    } //end if
  } //end if

  A11_->initializeProblemValues();
  A22_->setFieldOnDensityIsLinear(isLinear_);  // Set current state of linearity for F
  A22_->initializeProblemValues();

} //end initializeProblemValues
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, MatScalar value)
{

  LocalOrdinal schurBlockRow = physicsIdToSchurBlockId_[ownedPhysicsID];
  LocalOrdinal schurBlockCol = physicsIdToSchurBlockId_[boxPhysicsID];
  GlobalOrdinal rowGID = BLPM::ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = BLPM::boxToSolverGID(boxPhysicsID, boxNode);

  LocalOrdinal schurBlockNumber = 10 * schurBlockRow + schurBlockCol;

  switch (schurBlockNumber)
  {
  case 11:
    // A11 block
    A11_->insertMatrixValue(solverOrdering_[ownedPhysicsID], ownedMap_->getGlobalElement(ownedNode), rowGID, colGID, value);
    break;
  case 22:
    // A22 block
    // if poisson then blockColFlag = 0
    // if density then blockColFlag = 1
    // if cms then blockColFlag = 2
    if (isCmsEquation_[boxPhysicsID]) {
      A22_->insertMatrixValue(rowGID, colGID, value, 2);
    }else if (isDensityEquation_[boxPhysicsID]) {
      A22_->insertMatrixValue(rowGID, colGID, value, 1);
    }else if (isPoissonEquation_[boxPhysicsID]) {
      A22_->insertMatrixValue(rowGID, colGID, value, 0);
    }else{
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, "Unknown box physics ID in A22.");
    }
    break;
  case 21:
    // A12 block
    if (firstTime_) {
      if (rowGID!=curRowA21_) {
	// Dump the current contents of curRowValues_ into matrix and clear map
	insertRowA21();
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
	// Dump the current contents of curRowValues_ into matrix and clear map
	insertRowA12();
	curRowA12_=rowGID;
      }
      curRowValuesA12_[colGID] += value;
    }
    else
      A12Static_->sumIntoGlobalValues(rowGID, Array<GlobalOrdinal>(1,colGID), Array<MatScalar>(1,value));
    break;
  }

  if (debug_)
  {
    BLPM::insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
  }
} //insertMatrixValue
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
insertRowA12
()
{
  if (curRowValuesA12_.empty())
  {
    return;
  } //end if
  size_t numEntries = curRowValuesA12_.size();
  if (numEntries>indicesA12_.size())
  {
    indicesA12_.resize(numEntries);
    valuesA12_.resize(numEntries);
  } //end if
  LocalOrdinal i=0;

  for (ITER pos = curRowValuesA12_.begin(), e = curRowValuesA12_.end(); pos != e; ++pos)
  {
    indicesA12_[i] = pos->first;
    valuesA12_[i++] = pos->second;
  } //end for
   A12_->insertGlobalValues(curRowA12_, indicesA12_, valuesA12_);

  indicesA12_.clear();
  valuesA12_.clear();
  curRowValuesA12_.clear();

} //end insertRowA12
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
insertRowA21
()
{
  if (curRowValuesA21_.empty())
  {
    return;
  } //end if
  size_t numEntries = curRowValuesA21_.size();
  if (numEntries>indicesA21_.size())
  {
    indicesA21_.resize(numEntries);
    valuesA21_.resize(numEntries);
  } //end if
  LocalOrdinal i=0;

  for (ITER pos = curRowValuesA21_.begin(), e = curRowValuesA21_.end(); pos != e; ++pos)
  {
    indicesA21_[i] = pos->first;
    valuesA21_[i++] = pos->second;
  } //end for
   A21_->insertGlobalValues(curRowA21_, indicesA21_, valuesA21_);

  indicesA21_.clear();
  valuesA21_.clear();
  curRowValuesA21_.clear();
} //end insertRowA21
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
  {
    return; // nothing to do
  } //end if

  if (firstTime_)
  {
    insertRowA12(); // Dump any remaining entries
    insertRowA21(); // Dump any remaining entries

    RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
    pl->set( "Preserve Local Graph", true );
    if(!A12_->isFillComplete()){
      A12_->fillComplete(block2RowMap_,block1RowMap_,pl);
    }
    if(!A21_->isFillComplete()) {
      A21_->fillComplete(block1RowMap_,block2RowMap_,pl);
    }

    ArrayRCP<size_t> numEntriesPerRowA12(A12_->getRowMap()->getNodeNumElements());
    for (LocalOrdinal i = 0; i < A12_->getRowMap()->getNodeNumElements(); ++i) {
      numEntriesPerRowA12[i] = A12_->getNumEntriesInLocalRow( i );
    }
    ArrayRCP<size_t> numEntriesPerRowA21(A21_->getRowMap()->getNodeNumElements());
    for (LocalOrdinal i = 0; i < A21_->getRowMap()->getNodeNumElements(); ++i) {
      numEntriesPerRowA21[i] = A21_->getNumEntriesInLocalRow( i );
    }

    A12Graph_ = rcp(new GRAPH(A12_->getRowMap(), A12_->getColMap(), numEntriesPerRowA12, Tpetra::StaticProfile));
    A21Graph_ = rcp(new GRAPH(A21_->getRowMap(), A21_->getColMap(), numEntriesPerRowA21, Tpetra::StaticProfile));
    for (LocalOrdinal i = 0; i < A12_->getRowMap()->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      A12_->getLocalRowView( i, indices, values );
      A12Graph_->insertLocalIndices( i, indices );
    }
    for (LocalOrdinal i = 0; i < A21_->getRowMap()->getNodeNumElements(); ++i) {
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

    for (LocalOrdinal i = 0; i < A12_->getRowMap()->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      A12_->getLocalRowView( i, indices, values );
      A12Static_->sumIntoLocalValues( i, indices(), values() );
    }
    A12Static_->fillComplete(block2RowMap_,block1RowMap_,pl);
    for (LocalOrdinal i = 0; i < A21_->getRowMap()->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      A21_->getLocalRowView( i, indices, values );
      A21Static_->sumIntoLocalValues( i, indices(), values() );
    }
    A21Static_->fillComplete(block1RowMap_,block2RowMap_,pl);

    schurOperator_ = rcp(new ScTO(A11_, A12Static_, A21Static_, A22_));

    if (debug_)
    {
      globalMatrix_->fillComplete();
    } //end if

  } //end if
  RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
  if(!A12Static_->isFillComplete()){
    A12Static_->fillComplete(block2RowMap_,block1RowMap_,pl);
  }
  if(!A21Static_->isFillComplete()) {
    A21Static_->fillComplete(block1RowMap_,block2RowMap_,pl);
  }

  A11_->finalizeProblemValues();
  A22_->finalizeProblemValues();

  //  Check(true);
  isLinearProblemSet_ = true;
  firstTime_ = false;

} //end finalizeProblemValues
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
setupSolver
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error,
		     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");
#endif

  schurOperator_->ComputeRHS(*rhs1_, *rhs2_, *rhsSchur_);

#ifdef SUPPORTS_STRATIMIKOS
  thyraRhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(rhsSchur_);
  thyraLhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(lhs2_);
  thyraOp_ = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(schurOperator_);

  solver_ = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver_->readParameters(out.get());

  lowsFactory_ = solver_->createLinearSolveStrategy("");
  lows_ = linearOpWithSolve<Scalar>(*lowsFactory_, thyraOp_);
#else

  problem_ = rcp(new LinPROB(schurOperator_, lhs2_, rhsSchur_));
  int precond  = parameterList_->template get<int>( "Precond" );
  if (precond != AZ_none) {
    problem_->setLeftPrec(A22precond_);
  }

  TEUCHOS_TEST_FOR_EXCEPT(problem_->setProblem() == false);
  RCP<ParameterList> belosList = rcp(new ParameterList(parameterList_->sublist("belosList")));
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem_, belosList));
#endif

} //end setupSolver
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
solve
()
{
  //  A11_->Check(true);
#ifdef KDEBUG
  printf("\n\n\n\ndft_PolyLinProbMgr::solve()\n\n\n\n");
#endif

  // Solve the Schur complement system
#ifdef SUPPORTS_STRATIMIKOS
  SolveStatus<double> status = lows_->solve(Thyra::NOTRANS, *thyraRhs_, thyraLhs_.ptr());
#else
  try {
    ReturnType ret = solver_->solve();
  }
  catch (Belos::StatusTestError& e) {
    std::cout << "Belos failed to solve the linear problem! Belos threw exception "
	      << e.what() << std::endl;
  }
#endif

 // Compute the rest of the solution
  schurOperator_->ComputeX1(*rhs1_, *lhs2_, *lhs1_);

  if (debug_)
  {
    RCP<VEC> tmpRhs = rcp(new VEC(globalRowMap_));
    RCP<VEC> tmprhs1 = tmpRhs->offsetViewNonConst(block1RowMap_, 0)->getVectorNonConst(0);
    RCP<VEC> tmprhs2 = tmpRhs->offsetViewNonConst(block2RowMap_, block1RowMap_->getNodeNumElements())->getVectorNonConst(0);

    schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *tmprhs1, *tmprhs2);

    tmpRhs->update(-STS::one(), *globalRhs_, STS::one());
    Scalar resid = STS::zero();
    resid = tmpRhs->norm2();
    std::cout << "Global Residual for solution = " << resid << std::endl;
    bool writeMatrixNow = false;
    if (writeMatrixNow)
    {
      writeMatrix("A.dat", "GlobalMatrix", "GlobalMatrix");
      BLPM::writeLhs("x.dat");
      BLPM::writeRhs("b.dat");
      BLPM::writePermutation("p.dat");
    } //end if
  } //end if

} //end Solve
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
applyMatrix
(const ArrayView<const ArrayView<const Scalar> >& x) const
{
  this->setLhs(x);

  schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *rhs1_, *rhs2_);

  return (getRhs());
} //end applyMatrix
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
Check
(bool verbose) const
{
  A11_->Check(verbose);
  A22_->Check(verbose);
} //end Check
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>::
writeMatrix
(const char * filename, const char * matrixName, const char * matrixDescription) const
{
  if (debug_)
  {
    //return(EpetraExt::RowMatrixToMatrixMarketFile
    //      (filename, *globalMatrix_, matrixName, matrixDescription));
  } //end if
  else
  {
    return; // Not available if not in debug mode
  } //end else
} //end writeMatrix

TRAMONTO_INST_HELPER(dft_PolyLinProbMgr)
