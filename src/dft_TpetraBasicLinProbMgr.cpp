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

#include <MatrixMarket_Tpetra.hpp>

#include "dft_TpetraBasicLinProbMgr.hpp"

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_BasicLinProbMgr
(size_t numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm, RCP<Node> node)
  : isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    machineParamsSet_(false),
    groupByPhysics_(true),
    firstTime_(true),
    numUnknownsPerNode_(numUnknownsPerNode),
    numOwnedNodes_(0),
    numBoxNodes_(0),
    numGlobalNodes_(0),
    numGlobalBoxNodes_(0),
    numCoarsenedNodes_(0),
    numGlobalCoarsenedNodes_(0),
    comm_(comm),
    node_(node),
    curRow_(-1)
{
  // Setup global parameter list
  RCP<Tpetra::ParameterListConverter<Scalar,LocalOrdinal,GlobalOrdinal,Node> > pListConverter = rcp(new Tpetra::ParameterListConverter<Scalar,LocalOrdinal,GlobalOrdinal,Node>(parameterList));
  pListConverter->convert();
  ParameterList convertedParameters = pListConverter->getConvertedList();
  parameterList_ = rcp( new Teuchos::ParameterList(convertedParameters) );

  return;
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_BasicLinProbMgr()
{
  return;
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalRowMap
(const ArrayView<const GlobalOrdinal>& GIDs, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz)
{
  if (numGlobalNodes_!=0)
  {
    // Already been here
    return;
  }

  numOwnedNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
				  &numOwnedNodes_, &numGlobalNodes_);

  ownedMap_ = rcp(new MAP(numGlobalNodes_, GIDs, 0, comm_, node_));
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalColMap
(const ArrayView<const GlobalOrdinal> &GIDs, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz)
{
  if (numGlobalBoxNodes_!=0)
  {
    // Already been here
    return;
  }

  numBoxNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
				  &numBoxNodes_, &numGlobalBoxNodes_);

  boxMap_ = rcp(new MAP(numGlobalBoxNodes_, GIDs, 0, comm_, node_));
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
setCoarsenedNodesList(const ArrayView<const GlobalOrdinal> &GIDs)
{
  if (numGlobalCoarsenedNodes_!=0)
  {
    // Already been here
    return;
  }

  numCoarsenedNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
				  &numCoarsenedNodes_, &numGlobalCoarsenedNodes_);

  coarsenedNodesMap_ = rcp(new MAP(numGlobalCoarsenedNodes_, GIDs, 0, comm_, node_));
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeBlockStructure
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(isBlockStructureSet_, std::runtime_error, "Already set block structure.\n");
#endif

  const size_t numUnks = numOwnedNodes_*numUnknownsPerNode_;
  Array<LocalOrdinal> globalGIDList(numUnks);

  // Physics ordering for Basic Linear Problem is natural ordering:
  physicsOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); ++i)
  {
    physicsOrdering_[i] = i;
  }

  // Create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); ++i)
  {
    solverOrdering_[physicsOrdering_[i]]=i;
  }

  // Sanity check of physics ordering
  this->checkPhysicsOrdering();

  ArrayView<const GlobalOrdinal> GIDs = ownedMap_->getNodeElementList();
  LocalOrdinal k=0;
  if (groupByPhysics_)
  {
    for (LocalOrdinal i=0; i<numUnknownsPerNode_; ++i)
    {
      LocalOrdinal ii=physicsOrdering_[i];
      for (LocalOrdinal j=0; j<numOwnedNodes_; ++j)
      {
	globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
      }
    }
  }
  else
  {
    for (LocalOrdinal j=0; j<numOwnedNodes_; ++j)
    {
      for (LocalOrdinal i=0; i<numUnknownsPerNode_; ++i)
      {
	LocalOrdinal ii=physicsOrdering_[i];
	globalGIDList[k++] = ii + GIDs[j]*numUnknownsPerNode_;
      }
    }
  }

  Tpetra::global_size_t INVALID = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  globalRowMap_ = rcp(new MAP(INVALID, globalGIDList, 0, comm_, node_));
  globalMatrix_ = rcp(new MAT(globalRowMap_, 0));
  globalMatrix_->setObjectLabel("BasicLinProbMgr::globalMatrix");
  globalRhs_ = rcp(new VEC(globalRowMap_));
  globalLhs_ = rcp(new VEC(globalRowMap_));
  ownedToBoxImporter_ = rcp(new IMP(ownedMap_, boxMap_));

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;
}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(!isBlockStructureSet_, std::logic_error,
		     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
  TEUCHOS_TEST_FOR_EXCEPTION(!isGraphStructureSet_, std::logic_error,
		     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
#endif

  // Re-initialize the linear problem
  isLinearProblemSet_ = false;

  // We need to initialize the matrix the first time it is filled because matrix
  // entries may be put in on residual-only fills, which can occur before matrix
  // fills.

  if (firstTime_) {
    globalMatrix_->resumeFill();
    globalMatrix_->setAllToScalar(STMS::zero());
  } else {
    globalMatrixStatic_->resumeFill();
    globalMatrixStatic_->setAllToScalar(STMS::zero());
  }
  globalRhs_->putScalar(STS::zero());
  globalLhs_->putScalar(STS::zero());

 }

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRhsValue
(GlobalOrdinal ownedPhysicsID, GlobalOrdinal ownedNode, Scalar value)
{
  LocalOrdinal rhsLID = this->ownedToSolverLID(ownedPhysicsID, ownedNode); // Get solver LID
  globalRhs_->sumIntoLocalValue(rhsLID, value);
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, MatScalar value)
{
  GlobalOrdinal rowGID = this->ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = this->boxToSolverGID(boxPhysicsID, boxNode);

  if (firstTime_) {
    if (rowGID!=curRow_) {
      // Insert the current row values into the matrix and move on to the next row
      insertRow();
      curRow_=rowGID;
    }
    curRowValues_[colGID] += value;
  }
  else
    globalMatrixStatic_->sumIntoGlobalValues(rowGID, Array<GlobalOrdinal>(1,colGID), Array<MatScalar>(1,value));

}
//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  if (curRowValues_.empty()) return;

  for (ITER pos = curRowValues_.begin(), e = curRowValues_.end(); pos != e; ++pos) {
    indices_.append(pos->first);
    values_.append(pos->second);
  }

  globalMatrix_->insertGlobalValues(curRow_, indices_, values_);

  indices_.clear();
  values_.clear();
  curRowValues_.clear();
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
MatScalar
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
getMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode)
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(globalMatrixStatic_.get()==0, std::runtime_error, "Global Matrix is not constructed, must set debug flag to enable this feature.\n");
#endif

  GlobalOrdinal rowGID = this->ownedToSolverGID(ownedPhysicsID, ownedNode);
  GlobalOrdinal colGID = this->boxToSolverGID(boxPhysicsID, boxNode);
  size_t numEntries;
  ArrayView<const GlobalOrdinal> indices;
  ArrayView<const MatScalar> values;

  if (globalMatrixStatic_->isGloballyIndexed())
  {
    globalMatrixStatic_->getGlobalRowView(rowGID, indices, values);
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; ++i)
    {
      if (colGID==indices[i])
      {
	return(values[i]);
      }
    }
  }
  else
  {
    LocalOrdinal rowLID = globalMatrixStatic_->getRowMap()->getLocalElement(rowGID);
    LocalOrdinal colLID = globalMatrixStatic_->getColMap()->getLocalElement(colGID);
    globalMatrixStatic_->getLocalRowView(rowLID, indices, values);
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; ++i)
    {
      if (colLID==indices[i])
      {
	return(values[i]);
      }
    }
  }

  return(STMS::zero());
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID,
 const ArrayView<const LocalOrdinal>& boxNodeList, const ArrayView<const MatScalar>& values)
{
  for (LocalOrdinal i=0; i<boxNodeList.size(); ++i)
  {
    this->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNodeList[i], values[i]);
  }
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, const ArrayView<const LocalOrdinal> &boxPhysicsIDList,
 LocalOrdinal boxNode, const ArrayView<const MatScalar> &values)
{
  for (LocalOrdinal i=0; i<boxPhysicsIDList.size(); ++i)
  {
    this->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsIDList[i], boxNode, values[i]);
  }
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
    return;

  if (firstTime_) {

    insertRow();

    if(!globalMatrix_->isFillComplete()){
      RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
      pl->set( "Preserve Local Graph", true );
      globalMatrix_->fillComplete( pl );
    }
    ArrayRCP<size_t> numEntriesPerRow(globalRowMap_->getNodeNumElements());
    for (LocalOrdinal i = 0; i < globalRowMap_->getNodeNumElements(); ++i) {
      numEntriesPerRow[i] = globalMatrix_->getNumEntriesInLocalRow( i );
    }
    globalGraph_ = rcp(new GRAPH(globalRowMap_, globalMatrix_->getColMap(), numEntriesPerRow, Tpetra::StaticProfile));
    for (LocalOrdinal i = 0; i < globalRowMap_->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      globalMatrix_->getLocalRowView( i, indices, values );
      globalGraph_->insertLocalIndices( i, indices );
    }
    globalGraph_->fillComplete();
    globalMatrixStatic_ = rcp(new MAT(globalGraph_));
    globalMatrixStatic_->setAllToScalar(STMS::zero());
    for (LocalOrdinal i = 0; i < globalRowMap_->getNodeNumElements(); ++i) {
      ArrayView<const GlobalOrdinal> indices;
      ArrayView<const MatScalar> values;
      globalMatrix_->getLocalRowView( i, indices, values );
      globalMatrixStatic_->sumIntoLocalValues( i, indices(), values() );
    }
    globalMatrixStatic_->fillComplete();
    globalMatrixOp_ = rcp(new MMOP(globalMatrixStatic_));
  }

  if(!globalMatrixStatic_->isFillComplete()){
    RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
    globalMatrixStatic_->fillComplete( pl );
  }

  isLinearProblemSet_ = true;
  firstTime_ = false;

}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
setRhs
(const ArrayView<const ArrayView<const Scalar> >& b)
{
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; ++i) {
    for (LocalOrdinal j=0; j<numOwnedNodes_; ++j) {
      globalRhs_->replaceLocalValue(ownedToSolverLID(i,j), b[i][j]);
    }
  }
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
setLhs
(const ArrayView<const ArrayView<const Scalar> > &x) const
{
  ArrayRCP<const Scalar> xtmp;
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; ++i) {
    xtmp = exportC2R(Teuchos::arcpFromArrayView<const Scalar>(x[i]));
    for (LocalOrdinal j=0; j<numOwnedNodes_; ++j) {
      globalLhs_->replaceLocalValue(ownedToSolverLID(i,j), xtmp[j]);
    }
  }
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
getLhs
() const
{
  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalLhs_->get1dViewNonConst();

  for (LocalOrdinal i=0; i<numUnknownsPerNode_; ++i) {
    ArrayRCP<Scalar> temp(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; ++j) {
      temp[j] = tmp[ownedToSolverLID(i,j)];
    }
    ArrayOfPtrs[i] = this->importR2C(temp.getConst());
  }
  return ArrayOfPtrs;
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
getRhs
() const
{
  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalRhs_->get1dViewNonConst();
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; ++i) {
    ArrayOfPtrs[i] = Teuchos::arcp<Scalar>(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; ++j) {
      ArrayOfPtrs[i][j] = tmp[ownedToSolverLID(i,j)];
    }
  }
  return ArrayOfPtrs;
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
setMachineParams
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error,
		     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");
#endif

  // Get machine parameters
  n_ = globalMatrixStatic_->getGlobalNumCols();
  eps_ = Teuchos::ScalarTraits<Scalar>::eps();
  epsHalf_ = Teuchos::ScalarTraits<halfScalar>::eps();
  anorm_ = globalMatrixStatic_->getFrobeniusNorm();
  nae_ = Teuchos::as<Scalar>(n_) * anorm_ * eps_;
  snae_ = sqrt(Teuchos::as<Scalar>(n_))*anorm_*eps_;
  naeHalf_ = Teuchos::as<halfScalar>(Teuchos::as<Scalar>(n_)*anorm_*epsHalf_);
  snaeHalf_ = Teuchos::as<halfScalar>(sqrt(Teuchos::as<Scalar>(n_))*anorm_*epsHalf_);
  machineParamsSet_ = true;
  //  std::cout << "PROBLEM CONSTANTS--- anorm = " << anorm_ << " n = " << n_ << " snae = " << snae_ << std::endl;

}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
setupSolver
()
{
#ifdef KDEBUG
  TEUCHOS_TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error,
		     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");
#endif

  // Setup machine constants
  //  setMachineParams();

  scaling_ = parameterList_->template get<int>( "Scaling" );

  // Perform scaling
  if (scaling_ != AZ_none) {
    scalingMatrix_ = rcp(new SCALE(globalMatrixStatic_));
    rowScaleFactors_ = rcp(new VEC_M(globalRowMap_));
    rowScaleFactorsScalar_ = rcp(new VEC(globalRowMap_));

    scalingMatrix_->getRowScaleFactors( rowScaleFactors_, 1 );
    globalMatrixStatic_->leftScale( *rowScaleFactors_ );

#if MIXED_PREC == 1

    RCP<Tpetra::MultiVectorConverter<MatScalar,Scalar,LocalOrdinal,GlobalOrdinal,Node> > mvConverter;
    // Convert rowScaleFactors to scalar precision
    mvConverter->convert( *rowScaleFactors_, *rowScaleFactorsScalar_ );

    globalRhs_->elementWiseMultiply( STS::one(), *rowScaleFactorsScalar_, *globalRhs_, STS::zero() );

#elif MIXED_PREC == 0
    globalRhs_->elementWiseMultiply( STS::one(), *rowScaleFactors_, *globalRhs_, STS::zero() );
#endif
  }

#ifdef SUPPORTS_STRATIMIKOS
  thyraRhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalRhs_);
  thyraLhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalLhs_);
  thyraOp_ = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalMatrixStatic_);

  solver_ = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver_->readParameters(out.get());

  lowsFactory_ = solver_->createLinearSolveStrategy("");
  lows_ = linearOpWithSolve<Scalar>(*lowsFactory_, thyraOp_);
#else

  problem_ = rcp(new LinPROB(globalMatrixOp_, globalLhs_, globalRhs_));

  int precond  = parameterList_->template get<int>( "Precond" );
  if (precond != AZ_none) {
    LocalOrdinal overlapLevel = 0;
    preconditioner_ = rcp(new SCHWARZ(globalMatrixStatic_,overlapLevel));
    preconditionerOp_ = rcp(new SCHWARZ_OP(preconditioner_));
    ParameterList ifpack2List = parameterList_->sublist("ifpack2List");
    preconditioner_->setParameters(ifpack2List);
    preconditioner_->initialize();
    preconditioner_->compute();
    problem_->setLeftPrec(preconditionerOp_);
  }

  TEUCHOS_TEST_FOR_EXCEPT(problem_->setProblem() == false);
  RCP<ParameterList> belosList = rcp(new ParameterList(parameterList_->sublist("belosList")));
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem_, belosList));
#endif
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
solve
()
{
#ifdef KDEBUG
  printf("\n\n\n\ndft_BasicLinProbMgr::solve()\n\n\n\n");
#endif

  // Solve the global system
#ifdef SUPPORTS_STRATIMIKOS
  SolveStatus<Scalar> status = lows_->solve(Thyra::NOTRANS, *thyraRhs_, thyraLhs_.ptr());
#else
  try {
    solver_->solve();
  }
  catch (Belos::StatusTestError& e) {
    std::cout << "Belos failed to solve the linear problem! Belos threw exception "
	      << e.what() << std::endl;
  }
#endif

  bool writeMatrixNow = false;
  if (writeMatrixNow)
  {
    writeMatrix("A_ref.dat", "GlobalMatrix", "GlobalMatrix");
    writeLhs("x_ref.dat");
    writeRhs("b_ref.dat");
    writePermutation("p_ref.dat");
  }

  // Undo scaling
  if (scaling_ != AZ_none) {
    rowScaleFactors_->reciprocal( *rowScaleFactors_ );
    globalMatrixStatic_->leftScale( *rowScaleFactors_ );
#if MIXED_PREC == 1

    RCP<Tpetra::MultiVectorConverter<MatScalar,Scalar,LocalOrdinal,GlobalOrdinal,Node> > mvConverter;
    // Convert rowScaleFactors to scalar precision
    mvConverter->convert( *rowScaleFactors_, *rowScaleFactorsScalar_ );

    globalRhs_->elementWiseMultiply( STS::one(), *rowScaleFactorsScalar_, *globalRhs_, STS::zero() );

#elif MIXED_PREC == 0
    globalRhs_->elementWiseMultiply( STS::one(), *rowScaleFactors_, *globalRhs_, STS::zero() );
#endif
  }

}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
applyMatrix
(const ArrayView<const ArrayView<const Scalar> >& x) const
{
  this->setLhs(x);

  globalMatrixOp_->apply(*globalLhs_, *globalRhs_);

  return(this->getRhs());
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const ArrayRCP<const ArrayRCP<const Scalar> >& xOwned) const
{
  ArrayRCP<ArrayRCP<Scalar> > my_xBox = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; ++i) {
    my_xBox[i] = this->importR2C(xOwned[i]);
  }

  return(my_xBox);
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar>
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const ArrayRCP<const Scalar> &aOwned) const
{
  RCP<VEC> owned = rcp(new VEC(ownedMap_, aOwned()));
  RCP<VEC>  box = rcp(new VEC(boxMap_));

  box->doImport(*owned, *ownedToBoxImporter_, Tpetra::INSERT);

  return (box->get1dViewNonConst());
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar>
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
exportC2R
(const ArrayRCP<const Scalar>& aBox) const
{
  RCP<VEC> owned =  rcp(new VEC(ownedMap_));
  RCP<VEC> box = rcp(new VEC(boxMap_, aBox()));

  owned->doExport(*box, *ownedToBoxImporter_, Tpetra::INSERT);

  return(owned->get1dViewNonConst());
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
writeMatrix
(const char * filename, const char * matrixName, const char * matrixDescription) const  {

  std::string str_filename(filename);
  std::string str_matrixName(matrixName);
  std::string str_matrixDescription(matrixDescription);

  Tpetra::MatrixMarket::Writer<MAT>::writeSparseFile(str_filename,globalMatrixStatic_,str_matrixName,str_matrixDescription);

}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
writeLhs
(const char * filename) const  {

  std::string str_filename(filename);
  std::string str_matrixName("LHS");
  std::string str_matrixDescription("LHS");

  //  Tpetra::MatrixMarket::Writer<MAT>::writeDenseFile(str_filename,globalLhs_,str_matrixName,str_matrixDescription);
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
writeRhs
(const char * filename) const  {

  std::string str_filename(filename);
  std::string str_matrixName("RHS");
  std::string str_matrixDescription("RHS");

  //  Tpetra::MatrixMarket::Writer<MAT>::writeDenseFile(str_filename,globalRhs_,str_matrixName,str_matrixDescription);
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
writePermutation
(const char * filename) const  {
  //int dft_BasicLinProbMgr::writePermutation(const char * filename) const  {
  return;
    //(EpetraExt::BlockMapToMatrixMarketFile(filename, *globalRowMap_, " ", " ", false));
}

//=============================================================================
template <class Scalar, class MatScalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node>::
checkPhysicsOrdering()
const  {

  TEUCHOS_TEST_FOR_EXCEPTION(physicsOrdering_.size()==0, std::runtime_error, "No unknowns are registered with this problem manager.\n");

  size_t numUnks = physicsOrdering_.size();
  Array<Scalar> tmp(numUnks);
  for (LocalOrdinal i=0; i<numUnks; ++i)
  {
    LocalOrdinal curID = physicsOrdering_[i];
    TEUCHOS_TEST_FOR_EXCEPTION(curID <0, std::runtime_error, "Invalid unknown number " << curID << " is less than 0.\n");
    TEUCHOS_TEST_FOR_EXCEPTION(curID>=numUnks, std::runtime_error, "Invalid unknown number " << curID << " is greater than or equal to the number of unknowns (" << numUnks << ").\n");
      tmp[curID] = tmp[curID]+1;
      // Increment counter for this ID (at the end each ID should appear exactly one time).
  }

  for (LocalOrdinal i=0; i<numUnks; ++i)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(tmp[i]==0, std::runtime_error, "Unknown number " << i << " is not present and should be.\n");
    TEUCHOS_TEST_FOR_EXCEPTION(tmp[i]>1, std::runtime_error, "Unknown number " << i << " is present " << tmp[i] << " times and should be present only once.\n");
  }

}

TRAMONTO_INST_HELPER(dft_BasicLinProbMgr)
