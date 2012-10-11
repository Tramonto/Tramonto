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

#include "dft_HardSphereLinProbMgr.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "AztecOO.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "dft_HardSphereA11_Epetra_Operator.hpp"
#include "dft_HardSphereA22_Epetra_Operator.hpp"
#include "dft_A22Matrix_Epetra_Operator.hpp"
#include "dft_Schur_Epetra_Operator.hpp"


//=============================================================================
/*dft_HardSphereLinProbMgr::dft_HardSphereLinProbMgr(int numUnknownsPerNode, int * solverOptions, double * solverParams, MPI_Comm comm, bool formSchurMatrix, bool debug) 
  : dft_BasicLinProbMgr(numUnknownsPerNode, solverOptions, solverParams, comm),
    isA22Diagonal_(false),
    formSchurMatrix_(formSchurMatrix),
    debug_(debug),
    curRowA12_(-1),
    curRowA21_(-1) {
  //debug_=true;
  return;
  }*/
//=============================================================================
dft_HardSphereLinProbMgr::dft_HardSphereLinProbMgr(int numUnknownsPerNode, Teuchos::ParameterList * parameterList, MPI_Comm comm, bool formSchurMatrix, bool debug)
  : dft_BasicLinProbMgr(numUnknownsPerNode, parameterList, comm),
    isA22Diagonal_(false),
    formSchurMatrix_(formSchurMatrix),
    debug_(debug),
    curRowA12_(-1),
    curRowA21_(-1) {
  return;
}
//=============================================================================
dft_HardSphereLinProbMgr::~dft_HardSphereLinProbMgr() {
   return;
}
//=============================================================================
int dft_HardSphereLinProbMgr::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning

  // Create importer to map from owned nodes to box nodes

  ownedToBoxImporter_ = Teuchos::rcp(new Epetra_Import(*boxMap_, *ownedMap_));

  if (numGlobalNodes_==0 ||
      numGlobalBoxNodes_==0 ||
      indNonLocalEquations_.Length()==0 ||
      depNonLocalEquations_.Length()<0  ||
      densityEquations_.Length()==0) return(-1); // Error: One or more set methods not called
  
	TEUCHOS_TEST_FOR_EXCEPTION((numGlobalNodes_==0 ||
					   numGlobalBoxNodes_==0 ||
					   indNonLocalEquations_.Length()==0 ||
					   depNonLocalEquations_.Length()<0  ||
					   densityEquations_.Length()==0), std::logic_error, 
					   "One or more set methods not called.");
	//Not checking if poissonEquations_.Length()==0 because don't HAVE to have Poisson equations      
	//Not checking if gInvEquations_.Length()==0 because we don't have to have G inv equations
	
	// Fill physics ordering vector with the concatenated contents of the IDs for all physics types
  // Load Schur block mappings
  physicsOrdering_.Size(numUnknownsPerNode_);
  physicsIdToSchurBlockId_.Size(numUnknownsPerNode_);
  int * ptr = physicsOrdering_.Values();
  for (int i=0; i<indNonLocalEquations_.Length(); i++) {
    *ptr++ = indNonLocalEquations_[i];
    physicsIdToSchurBlockId_[indNonLocalEquations_[i]] = 1;
  }
  for (int i=0; i<depNonLocalEquations_.Length(); i++) {
    *ptr++ = depNonLocalEquations_[i];
    physicsIdToSchurBlockId_[depNonLocalEquations_[i]] = 1;
  }
  for (int i=0; i<densityEquations_.Length(); i++) {
    *ptr++ = densityEquations_[i];
    physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
  }

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.Size(numUnknownsPerNode_);
  for (int i=0; i<physicsOrdering_.Length(); i++) solverOrdering_[physicsOrdering_[i]]=i;

  // Special setup for coarsened nodes
  // Build int vector of with 1 if node is coarsened
  ownedNodeIsCoarsened_ = Teuchos::rcp(new Epetra_IntVector(*ownedMap_)); // Assume no nodes are coarsened
  boxNodeIsCoarsened_ = Teuchos::rcp(new Epetra_IntVector(*boxMap_)); // Assume no nodes are coarsened
  
  if (numGlobalCoarsenedNodes_>0) {
    setA22BlockIsDiagonal(false); // A22 block is not diagonal when some nodes are coarsened
    for (int i=0; i<numCoarsenedNodes_; i++)
      (*ownedNodeIsCoarsened_)[ownedMap_->LID(coarsenedNodesMap_->GID(i))] = 1;
    boxNodeIsCoarsened_->Import(*ownedNodeIsCoarsened_, *ownedToBoxImporter_, Insert); // Now each processor knows which of its box nodes is coarsened
  }
  const int numUnks = numOwnedNodes_*numUnknownsPerNode_;
  const int numRealNodes = numOwnedNodes_ - numCoarsenedNodes_;
  const int numUnks1 = numRealNodes*(indNonLocalEquations_.Length()+depNonLocalEquations_.Length());
  const int numUnks2 = numRealNodes*(densityEquations_.Length()) + numCoarsenedNodes_ * numUnknownsPerNode_;
  assert(numUnks==(numUnks1+numUnks2));  // Sanity test
  const int numIndNonLocal = numRealNodes*(indNonLocalEquations_.Length());
  const int numDepNonLocal = numRealNodes*(depNonLocalEquations_.Length());
  const int numDensity = numRealNodes*(densityEquations_.Length());
  Epetra_IntSerialDenseVector globalGIDList(numUnks);

  int * GIDs = ownedMap_->MyGlobalElements();
  int k=0;
  int k1 = (numOwnedNodes_ - numCoarsenedNodes_) * numUnknownsPerNode_; // starting point for coarsened variables
  for (int i=0; i<numUnknownsPerNode_; i++) {
    int ii=physicsOrdering_[i];
    if (numCoarsenedNodes_==0) {
      for (int j=0; j<numOwnedNodes_; j++) 
	globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
    }
    else {
      for (int j=0; j<numOwnedNodes_; j++) {
	int curGID = GIDs[j];
	if (coarsenedNodesMap_->MyGID(curGID))
	  globalGIDList[k1++] = ii*numGlobalNodes_ + GIDs[j];
	else
	  globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
      }
    }
  }

  ptr = globalGIDList.Values();
  globalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks, ptr, 0, comm_));
  block1RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks1, ptr, 0, comm_));
  block2RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks2, ptr+numUnks1, 0, comm_));
  indNonLocalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numIndNonLocal, ptr, 0, comm_));
  depNonLocalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numDepNonLocal, ptr+numIndNonLocal, 0, comm_));
  /*
    std::cout << " Global Row Map" << *globalRowMap_ << std::endl
    << " Block 1     Row Map " << *block1RowMap_ << std::endl
    << " Block 2     Row Map " << *block2RowMap_ << std::endl
    << " DepNonLocal Row Map " << *depNonLocalRowMap_ << std::endl;
  */
  A11_ = Teuchos::rcp(new dft_HardSphereA11_Epetra_Operator(*indNonLocalRowMap_, *depNonLocalRowMap_, *block1RowMap_));
  A12_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *block1RowMap_, 0)); A12_->SetLabel("HardSphere::A12");
  A21_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *block2RowMap_, 0)); A21_->SetLabel("HardSphere::A21");
  if (isA22Diagonal_) {
    A22Diagonal_ = Teuchos::rcp(new dft_HardSphereA22_Epetra_Operator(*block2RowMap_));
  }
  else {
    A22Matrix_ = Teuchos::rcp(new dft_A22Matrix_Epetra_Operator(*block2RowMap_));
  }
  if (debug_) {
    globalMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalRowMap_, 0));
    globalMatrix_->SetLabel("HardSphere::globalMatrix");
    }
  else
    globalMatrix_ = Teuchos::null; // not used by this solver
  
  globalRhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalLhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));

  rhs1_ = Teuchos::rcp(new Epetra_Vector(View, *block1RowMap_, globalRhs_->Values()));
  rhs2_ = Teuchos::rcp(new Epetra_Vector(View, *block2RowMap_, globalRhs_->Values()+numUnks1));
  rhsSchur_ = Teuchos::rcp(new Epetra_Vector(*rhs2_));
  lhs1_ = Teuchos::rcp(new Epetra_Vector(View, *block1RowMap_, globalLhs_->Values()));
  lhs2_ = Teuchos::rcp(new Epetra_Vector(View, *block2RowMap_, globalLhs_->Values()+numUnks1));

  if (isA22Diagonal_) 
    schurOperator_ = Teuchos::rcp(new dft_Schur_Epetra_Operator(A11_.get(), A12_.get(), A21_.get(), A22Diagonal_.get()));
  else
    schurOperator_ = Teuchos::rcp(new dft_Schur_Epetra_Operator(A11_.get(), A12_.get(), A21_.get(), A22Matrix_.get()));

  implicitProblem_ = Teuchos::rcp(new Epetra_LinearProblem(schurOperator_.get(), lhs2_.get(), rhsSchur_.get()));

    

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::initializeProblemValues() {
  
	TEUCHOS_TEST_FOR_EXCEPTION(!isBlockStructureSet_, std::logic_error, 
					   "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
	TEUCHOS_TEST_FOR_EXCEPTION(!isGraphStructureSet_, std::logic_error, 
					   "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
	isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    A12_->PutScalar(0.0);
    A21_->PutScalar(0.0);
    globalRhs_->PutScalar(0.0);
    globalLhs_->PutScalar(0.0);
    if (debug_) globalMatrix_->PutScalar(0.0);
  }

  A11_->initializeProblemValues();
  if (isA22Diagonal_)
    A22Diagonal_->initializeProblemValues();
  else
    A22Matrix_->initializeProblemValues();
  
  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value) {

  bool schurBlockRow1 = (physicsIdToSchurBlockId_[ownedPhysicsID]==1 && (*ownedNodeIsCoarsened_)[ownedNode]==0);
  bool schurBlockCol1 = (physicsIdToSchurBlockId_[boxPhysicsID]==1 && (*boxNodeIsCoarsened_)[boxNode]==0);
  int rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  int colGID = boxToSolverGID(boxPhysicsID, boxNode);
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
      A21_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
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
      A12_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
  }

  if (debug_) dft_BasicLinProbMgr::insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
  
  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::insertRowA12() {
  if (curRowValuesA12_.empty()) return(0);
  int numEntries = curRowValuesA12_.size();
  if (numEntries>indicesA12_.Length()) {
    indicesA12_.Resize(numEntries);
    valuesA12_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos = curRowValuesA12_.begin(); pos != curRowValuesA12_.end(); ++pos) {
    indicesA12_[i] = pos->first;
    valuesA12_[i++] = pos->second;
  }
   A12_->InsertGlobalValues(curRowA12_, numEntries, valuesA12_.Values(), indicesA12_.Values());

  curRowValuesA12_.clear();
  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::insertRowA21() {
  if (curRowValuesA21_.empty()) return(0);
  int numEntries = curRowValuesA21_.size();
  if (numEntries>indicesA21_.Length()) {
    indicesA21_.Resize(numEntries);
    valuesA21_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos = curRowValuesA21_.begin(); pos != curRowValuesA21_.end(); ++pos) {
    indicesA21_[i] = pos->first;
    valuesA21_[i++] = pos->second;
  }
   A21_->InsertGlobalValues(curRowA21_, numEntries, valuesA21_.Values(), indicesA21_.Values());

  curRowValuesA21_.clear();
  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) {
    insertRowA12(); // Dump any remaining entries
    A12_->FillComplete(*block2RowMap_,*block1RowMap_);
    A12_->OptimizeStorage();
    insertRowA21(); // Dump any remaining entries
    A21_->FillComplete(*block1RowMap_,*block2RowMap_);
    A21_->OptimizeStorage();

    if (debug_) globalMatrix_->FillComplete();
  }
  //std::cout << *A12_ << endl 
  //          << *A21_ << endl;

  A11_->finalizeProblemValues();
  if (isA22Diagonal_)
    A22Diagonal_->finalizeProblemValues();
  else
    A22Matrix_->finalizeProblemValues();

  //Check(true);
  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::setupSolver() {

	
	TEUCHOS_TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error, 
					   "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");
	
	
  schurOperator_->ComputeRHS(*rhs1_, *rhs2_, *rhsSchur_);
  if (solver_ != Teuchos::null ) return(0);  //Already setup
	
  if (formSchurMatrix_) {// We have S explicitly available, so let's use it
    if (isA22Diagonal_)
      schurOperator_->SetSchurComponents(A11_->getA11invMatrix(), A22Diagonal_->getA22Matrix());
    else
      schurOperator_->SetSchurComponents(A11_->getA11invMatrix(), A22Matrix_->getA22Matrix());
    implicitProblem_->SetOperator(schurOperator_->getSchurComplement());
  }
  solver_ = Teuchos::rcp(new AztecOO(*implicitProblem_));

  solver_->SetParameters(*parameterList_);
  //  if (solverOptions_!=0) solver_->SetAllAztecOptions(solverOptions_);
  //  if (solverParams_!=0) solver_->SetAllAztecParams(solverParams_);

  //  const int * options = solver_->GetAllAztecOptions();
  //  const double * params = solver_->GetAllAztecParams();

  solver_->SetAztecOption(AZ_scaling, AZ_none); 
  int maxiter = 500;
  solver_->SetAztecOption(AZ_max_iter, maxiter);
  solver_->SetAztecOption(AZ_kspace, maxiter); 
  //solver_->SetAztecParam(AZ_tol, 1.0e-12); 
  //solver_->SetAztecOption(AZ_conv, AZ_noscaled); 
  if (isA22Diagonal_)
    solver_->SetPrecOperator(A22Diagonal_.get());
  else
    solver_->SetPrecOperator(A22Matrix_.get());
  //solver_->SetAztecParam(AZ_ill_cond_thresh, 0.0);
  //solver_->SetAztecOption(AZ_precond, AZ_none);
  //solver_->SetAztecOption(AZ_solver, AZ_bicgstab);
  //solver_->SetAztecOption(AZ_precond, AZ_dom_decomp);
  //solver_->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  //solver_->SetAztecParam(AZ_ilut_fill, 4.0);
  //solver_->SetAztecParam(AZ_drop, 0.0);

  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::solve() {
  
  //writeMatrix("2D.mm", "Small HardSpheremer Matrix", "Global Matrix from Small HardSpheremer Problem");
  //abort();

  //  const int * options = solver_->GetAllAztecOptions();
  //  const double * params = solver_->GetAllAztecParams();

  // solver_->Iterate(options[AZ_max_iter], params[AZ_tol]); // Try to solve
  solver_->Iterate(Teuchos::getParameter<int>(*parameterList_, "Max_iter"), Teuchos::getParameter<double>(*parameterList_, "Tol")); // Try to solve
  schurOperator_->ComputeX1(*rhs1_, *lhs2_, *lhs1_); // Compute rest of solution
  if (debug_) {
    Epetra_Vector tmpRhs(*globalRowMap_);
    Epetra_Vector tmprhs1(View, *block1RowMap_, tmpRhs.Values());
    Epetra_Vector tmprhs2(View, *block2RowMap_, tmpRhs.Values()+block1RowMap_->NumMyElements());
    
    schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, tmprhs1, tmprhs2);
    
    tmpRhs.Update(-1.0, *globalRhs_, 1.0);
    double resid=0.0;
    tmpRhs.Norm2(&resid);
    std::cout << "Global Residual for solution = " << resid << std::endl;
    bool writeMatrixNow = true;
    if (writeMatrixNow) {
      writeMatrix("A.dat", "GlobalMatrix", "GlobalMatrix");
      writeLhs("x.dat");
      writeRhs("b.dat");
      writePermutation("p.dat");
      abort();
    }
  }
  //std::cout << "Global RHS = " << *globalRhs_ << std::endl
  //          << "Global LHS = " << *globalLhs_ << std::endl;

  //solver_->AdaptiveIterate(solverOptions_[AZ_max_iter], 5, solverParams_[AZ_tol]); // Try to solve
  return(0);
}
//=============================================================================
int dft_HardSphereLinProbMgr::applyMatrix(const double** x, double** b) const {
  
  setLhs(x);
  schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *rhs1_, *rhs2_);
  getRhs(b);
  
  return(0);
}
//=============================================================================
  int dft_HardSphereLinProbMgr::Check(bool verbose) const {

  int ierr1 = A11_->Check(verbose);
  int ierr2;
  if (isA22Diagonal_)
    ierr2 = A22Diagonal_->Check(verbose);
  else
    ierr2 = A22Diagonal_->Check(verbose);
  if (ierr1!=0 || ierr2!=0) return(-1);
  return(0);
    
  }
//=============================================================================
int dft_HardSphereLinProbMgr::writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const  {
//  if (debug_){
    return(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_, matrixName, matrixDescription));
//  }
//  else{
 //   return(-1); // Not available if not in debug mode
  //}
}
