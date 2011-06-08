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

#include "dft_TpetraBasicLinProbMgr.hpp"

//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_BasicLinProbMgr
(size_t numUnknownsPerNode, RCP<ParameterList> parameterList, 
 RCP<const COMM> comm)
  : numUnknownsPerNode_(numUnknownsPerNode),
    parameterList_(parameterList),
    numOwnedNodes_(0),
    numBoxNodes_(0),
    numGlobalNodes_(0),
    numGlobalBoxNodes_(0),
    numCoarsenedNodes_(0),
    numGlobalCoarsenedNodes_(0),
    comm_(comm),
    isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    groupByPhysics_(true),
    firstTime_(true),
    curRow_(-1) 
{ 
#ifdef KDEBUG
  printf("\n\n\nCreated a BasicLinProbMgr.\n\n\n");
#endif
  return;
} //end constructor
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_BasicLinProbMgr
() 
{
  return;
} //end destructor
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalRowMap
(const ArrayView<const GlobalOrdinal>& GIDs, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz) 
{
  if (numGlobalNodes_!=0) 
  {
    return; // Already been here
  }
  numOwnedNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1, 
                               &numOwnedNodes_, &numGlobalNodes_);

  ownedMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs, 0, comm_));
} //end setNodalRowMap
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalColMap
(const ArrayView<const GlobalOrdinal> &GIDs, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz) 
{
  
  if (numGlobalBoxNodes_!=0) 
  {
    return; // Already been here
  } //end if

  numBoxNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1, 
                               &numBoxNodes_, &numGlobalBoxNodes_);

  boxMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs, 0, comm_));
  //std::cout << " Box Map" << *boxMap_ << std::endl;
} //end setNodalColMap
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setCoarsenedNodesList
(const ArrayView<const GlobalOrdinal> &GIDs) 
{
  
  if (numGlobalCoarsenedNodes_!=0) 
  {
    return; // Already been here
  } //end if

  numCoarsenedNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
                               &numCoarsenedNodes_, &numGlobalCoarsenedNodes_);

  coarsenedNodesMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs, 0, comm_));
  //std::cout << " Coarsened Nodes Map" << *coarsenedNodesMap_ << std::endl;
} //end setCarsenedNodesList
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeBlockStructure
() 
{
  TEST_FOR_EXCEPTION(isBlockStructureSet_, std::runtime_error, "Already set block structure.\n");
  
  const size_t numUnks = numOwnedNodes_*numUnknownsPerNode_;
  Array<LocalOrdinal> globalGIDList(numUnks);

  // Physics ordering for Basic Linear Problem is natural ordering:
  physicsOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++) 
  {
    physicsOrdering_[i] = i;
  } //end for

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++)
  {
    solverOrdering_[physicsOrdering_[i]]=i;
  } //end for

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  ArrayView<const GlobalOrdinal> GIDs = ownedMap_->getNodeElementList();
  LocalOrdinal k=0;
  if (groupByPhysics_) 
  {
    for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) 
    {
      LocalOrdinal ii=physicsOrdering_[i];
      for (LocalOrdinal j=0; j<numOwnedNodes_; j++) 
      {
	      globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
      } //end for
    } //end for
  } //end if
  else
  {
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) 
      {
        LocalOrdinal ii=physicsOrdering_[i];
        globalGIDList[k++] = ii + GIDs[j]*numUnknownsPerNode_;
      } //end for
    } //end for
  } //end else
  
  globalRowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList, 0, comm_));

  //std::cout << " Global Row Map" << *globalRowMap_ << std::endl;

  globalMatrix_ = rcp(new MAT(globalRowMap_, 0));
  globalMatrix_->setObjectLabel("BasicLinProbMgr::globalMatrix");
  globalRhs_ = rcp(new VEC(globalRowMap_));
  globalLhs_ = rcp(new VEC(globalRowMap_));
  //problem_ = rcp(new PROB(globalMatrix_, globalRhs_, globalLhs_));    

  ownedToBoxImporter_ = rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(ownedMap_, boxMap_));

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;
} //end finalizeBlockStructure
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
() 
{
  
  TEST_FOR_EXCEPTION(!isBlockStructureSet_, std::logic_error, 
                     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
  TEST_FOR_EXCEPTION(!isGraphStructureSet_, std::logic_error, 
                     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
	
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

 // AGS: I found that we needed to initialize the matrix even the 
 // first time a matrix was filled, because matrix entries are being put
 // in on residual-only fills, which can occur before matrix fills.

 if (!firstTime_) 
 {
    globalMatrix_->resumeFill();
    globalMatrix_->setAllToScalar(0.0);
    globalMatrix_->fillComplete();
    //    globalRhs_->putScalar(0.0);
    globalRhs_->putScalar(0.0);
    globalLhs_->putScalar(0.0);
  }
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRhsValue
(GlobalOrdinal ownedPhysicsID, GlobalOrdinal ownedNode, Scalar value) 
{

  LocalOrdinal rhsLID = ownedToSolverLID(ownedPhysicsID, ownedNode); // Get solver LID
  globalRhs_->sumIntoLocalValue(rhsLID, value);
  //cout << std::setprecision(2);
  //cout << "b[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode<<"] = " << value << endl;
} //end insertRhsValue
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, Scalar value) 
{

  GlobalOrdinal rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = boxToSolverGID(boxPhysicsID, boxNode);

  //cout << std::setprecision(2);
  //cout << "A[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode
  //     << "][boxPhysicsID="  <<boxPhysicsID  <<"][boxNode="  <<boxNode
  //     << "][rowGID="        <<rowGID        <<"][colGID="   <<colGID  
  //     << "] = " << value << endl;
  if (firstTime_) 
  {
    if (rowGID!=curRow_) 
    { 
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    } //end if
    curRowValues_[colGID] += value;
  } //end if
  else
  {
    Array<Scalar> vals(1);
    vals[0] = value;
    //    if(globalMatrix_->isGloballyIndexed()){
    Array<GlobalOrdinal> globalCols(1);
      globalCols[0] = colGID;
      //printf("Inserting global col[%d]\n", globalCols[0]);
      globalMatrix_->sumIntoGlobalValues(rowGID, globalCols, vals);
      //    } else {
      //      Array<LocalOrdinal> localCols(1);
      //      localCols[0] = globalMatrix_->getColMap()->getLocalElement(colGID);
      //      LocalOrdinal rowLID = globalMatrix_->getRowMap()->getLocalElement(rowGID);
      //printf("Inserting local col[%d]\n", localCols[0]);
      //      globalMatrix_->sumIntoLocalValues(rowLID, localCols, vals);
      //    }
  } //end else
} //end insertMatrixValue
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
() 
{
  if (curRowValues_.empty())
  {
    return;
  } //end if
  typename std::map<LocalOrdinal, Scalar>::iterator pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos) 
  {
    indices_.append(pos->first);
    values_.append(pos->second);
  } //end for
  //  if(globalMatrix_->isGloballyIndexed()){
    globalMatrix_->insertGlobalValues(curRow_, indices_, values_);
    //  } else {
    //    LocalOrdinal localRow = globalMatrix_->getRowMap()->getLocalElement(curRow_);
    //    Array<LocalOrdinal> localCol;
    //    localCol.resize(indices_.size());
    //    for(int i = 0; i < indices_.size(); i++){
    //      localCol[i] = globalMatrix_->getRowMap()->getLocalElement(indices_[i]);
      //TODO: Figure out why getColMap()->getLocalElement returns -1
    //    }
    //    globalMatrix_->resumeFill();
    //    globalMatrix_->insertLocalValues(localRow, localCol, values_);
    //    globalMatrix_->fillComplete();
    //  }

  indices_.clear();
  values_.clear();
  curRowValues_.clear();
} //end insertRow
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar 
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode) 
{

  TEST_FOR_EXCEPTION(globalMatrix_.get()==0, std::runtime_error, "Global Matrix is not constructed, must set debug flag to enable this feature.\n");
  
  GlobalOrdinal rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = boxToSolverGID(boxPhysicsID, boxNode);
  size_t numEntries;
  ArrayView<const GlobalOrdinal> indices;
  ArrayView<const Scalar> values;
  if (globalMatrix_->isGloballyIndexed()) 
  {
    globalMatrix_->getGlobalRowView(rowGID, indices, values); // get view of current row
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; i++)
    {
      if (colGID==indices[i]) 
      {
        return(values[i]);
      } //end if
    } //end for
  } //end if
  else 
  {
    rowGID = globalMatrix_->getRowMap()->getLocalElement(rowGID); // get local row ID
    colGID = globalMatrix_->getColMap()->getLocalElement(colGID); // get local column ID
    globalMatrix_->getLocalRowView(rowGID, indices, values); // get view of current row
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; i++)
    { 
      if (colGID==indices[i])
      {
        return(values[i]);
      } //end if
    } //end for
  } //end else
  
  return(0.0);
} //end getMatrixValue
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, const ArrayView<const LocalOrdinal>& boxNodeList, const ArrayView<const Scalar>& values) 
{
  
  for (LocalOrdinal i=0; i<boxNodeList.size(); i++)
  { 
    insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNodeList[i], values[i]);
  } //end for
} //end insertMatrixValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, const ArrayView<const LocalOrdinal> &boxPhysicsIDList, LocalOrdinal boxNode, const ArrayView<const Scalar> &values) 
{
  
  for (LocalOrdinal i=0; i<boxPhysicsIDList.size(); i++) 
  {
    insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsIDList[i], boxNode, values[i]);
  } //end for
} //end insertMatrixValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
() 
{
  if (isLinearProblemSet_)
  {
    return; // nothing to do
  } //end if

  if (firstTime_) 
  {
    insertRow();
    if(!globalMatrix_->isFillComplete()){
      globalMatrix_->fillComplete();
    }
  } //end if
  if(!globalMatrix_->isFillComplete()){
    globalMatrix_->fillComplete( );
  }

  isLinearProblemSet_ = true;
  firstTime_ = false;

  //writeMatrix("basica.dat", "", "");
} //end finalizeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setRhs
(const ArrayView<const ArrayView<const Scalar> >& b) 
{
  //printf("SETRHS called!\n");
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++)
  {
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      globalRhs_->replaceLocalValue(ownedToSolverLID(i,j), b[i][j]);
    } //end for
  } //end for
} //end setRHS
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setLhs
(const ArrayView<const ArrayView<const Scalar> > &x) const 
{
  //printf("calling setlhs\n");
  //Array<Scalar> xtmp(numOwnedNodes_); // Temp vector to hold local x values
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) 
  {
    ArrayRCP<const Scalar> xtmp = exportC2R(Teuchos::arcpFromArrayView<const Scalar>(x[i])); // Use simple import
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      globalLhs_->replaceLocalValue(ownedToSolverLID(i,j), xtmp[j]);
    } //end for
  } //end for
} //end setLHS
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getLhs
() const 
{

  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalLhs_->get1dViewNonConst();
  
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) 
  {
    ArrayRCP<Scalar> temp(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      temp[j] = tmp[ownedToSolverLID(i,j)];
    } //end for
    ArrayOfPtrs[i] = importR2C(temp.getConst()); // Use simple import (uses view so doesn't work)
  } //end for
  return ArrayOfPtrs;
} //end getLHS
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getRhs
() const 
{

  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalRhs_->get1dViewNonConst();
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++)
  {
    ArrayOfPtrs[i] = Teuchos::arcp<Scalar>(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      ArrayOfPtrs[i][j] = tmp[ownedToSolverLID(i,j)];
    } //end for
  } //end for
  return ArrayOfPtrs;
} //end getRHS
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setupSolver
() 
{
  TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error, 
                     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");

#ifdef SUPPORTS_STRATIMIKOS
  thyraRhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalRhs_);
  thyraLhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalLhs_);
  thyraOp_ = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalMatrix_);

  solver_ = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver_->readParameters(out.get());

  lowsFactory_ = solver_->createLinearSolveStrategy("");
  lows_ = linearOpWithSolve<Scalar>(*lowsFactory_, thyraOp_);
#else
  //printf("Node num rows = %d\n", globalMatrix_->getNodeNumRows());
  //printf("Node num cols = %d\n", globalMatrix_->getNodeNumCols());
  problem_ = rcp(new LinPROB(globalMatrix_, globalLhs_, globalRhs_));
  preconditioner_ = Ifpack2::Factory::create<MAT>("ILUT", globalMatrix_);
  preconditioner_->setParameters(*parameterList_);
  preconditioner_->initialize();
  preconditioner_->compute();
  problem_->setLeftPrec(preconditioner_);
  TEST_FOR_EXCEPT(problem_->setProblem() == false);
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem_, parameterList_));
#endif
} //end setupSolver
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
solve
() 
{
#ifdef KDEBUG
  printf("\n\n\n\ndft_BasicLinProbMgr::solve()\n\n\n\n");
#endif

#ifdef SUPPORTS_STRATIMIKOS
  SolveStatus<Scalar> status = lows_->solve(Thyra::NOTRANS, *thyraRhs_, thyraLhs_.ptr());
#else
  ReturnType ret = solver_->solve();
#endif

  bool writeMatrixNow = false;
  if (writeMatrixNow) 
  {
    writeMatrix("A_ref.dat", "GlobalMatrix", "GlobalMatrix");
    writeLhs("x_ref.dat");
    writeRhs("b_ref.dat");
    writePermutation("p_ref.dat");
  }
} //end solve
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyMatrix
(const ArrayView<const ArrayView<const Scalar> >& x) const 
{
  setLhs(x);
  //printf("Apply called!\n");
  globalMatrix_->apply(*globalLhs_, *globalRhs_);
  return(getRhs());
} //end applyMatrix
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const ArrayRCP<const ArrayRCP<const Scalar> >& xOwned) const 
{
  ArrayRCP<ArrayRCP<Scalar> > xBox = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) 
  {
    xBox[i] = importR2C(xOwned[i]);
  } //end for
  return(xBox);
} //end importR2C
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar> 
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const ArrayRCP<const Scalar> &aOwned) const 
{
  
  RCP<VEC> owned = rcp(new VEC(ownedMap_, aOwned()));
  RCP<VEC>  box = rcp(new VEC(boxMap_));
  
  
  //std::cout << " TargetMap " << ownedToBoxImporter_->getTargetMap() << std::endl;
  //std::cout << " SourceMap " << ownedToBoxImporter_->getSourceMap() << std::endl;
  //std::cout << " boxMap " << boxMap_ << std::endl;
  //std::cout << " ownedMap " << ownedMap_ << std::endl;

  box->doImport(*owned, *ownedToBoxImporter_, Tpetra::INSERT);

  return (box->get1dViewNonConst());
} //end importR2C
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar> 
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
exportC2R
(const ArrayRCP<const Scalar>& aBox) const 
{
  
  RCP<VEC> owned =  rcp(new VEC(ownedMap_));
  RCP<VEC> box = rcp(new VEC(boxMap_, aBox()));
  
  owned->doExport(*box, *ownedToBoxImporter_, Tpetra::INSERT); // Use importer, but zero out off-processor contributions.
  return(owned->get1dViewNonConst());
} //end exportC2R
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeMatrix
(const char * filename, const char * matrixName, const char * matrixDescription) const  
{
  return;//(EpetraExt::RowMatrixToMatrixMarketFile(filename, 
         //*globalMatrix_, matrixName, matrixDescription));
} //end writeMatrix
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeLhs
(const char * filename) const  
{
  return;//(EpetraExt::MultiVectorToMatlabFile(filename, *globalLhs_));
} //end writeLhs
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeRhs
(const char * filename) const  
{
  return;//(EpetraExt::MultiVectorToMatlabFile(filename, *globalRhs_));
} //end writeRhs
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writePermutation
(const char * filename) const  
{
  return;//(EpetraExt::BlockMapToMatrixMarketFile(filename, *globalRowMap_, " ", " ", false));
} //end writePermuation
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
checkPhysicsOrdering() const  
{
  TEST_FOR_EXCEPTION(physicsOrdering_.size()==0, std::runtime_error, "No unknowns are registered with this problem manager.\n"); 

  size_t numUnks = physicsOrdering_.size();
  Array<Scalar> tmp(numUnks);
  for (LocalOrdinal i=0; i<numUnks; i++) 
  {
    LocalOrdinal curID = physicsOrdering_[i];
    TEST_FOR_EXCEPTION(curID <0, std::runtime_error, "Invalid unknown number " << curID << " is less than 0.\n");
    TEST_FOR_EXCEPTION(curID>=numUnks, std::runtime_error, "Invalid unknown number " << curID << " is greater than or equal to the number of unknowns (" << numUnks << ").\n");
      tmp[curID] = tmp[curID]+1; 
      // Increment counter for this ID (at the end each ID should appear exactly one time).
  } //end for

  for (LocalOrdinal i=0; i<numUnks; i++) 
  {
    TEST_FOR_EXCEPTION(tmp[i]==0, std::runtime_error, "Unknown number " << i << " is not present and should be.\n");
    TEST_FOR_EXCEPTION(tmp[i]>1, std::runtime_error, "Unknown number " << i << " is present " << tmp[i] << " times and should be present only once.\n");
  } //end for
} //end checkPhysicsOrdering
template class dft_BasicLinProbMgr<double, int, int>;
