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

#include "dft_PolyLinProbMgr.hpp"
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
#include "dft_PolyA11_Epetra_Operator.hpp"
#include "dft_PolyA11_Coulomb_Epetra_Operator.hpp"
#include "dft_PolyA22_Epetra_Operator.hpp"
#include "dft_PolyA22_Coulomb_Epetra_Operator.hpp"
#include "dft_Schur_Epetra_Operator.hpp"


//=============================================================================
/*dft_PolyLinProbMgr::dft_PolyLinProbMgr(int numUnknownsPerNode, int * solverOptions, double * solverParams, MPI_Comm comm, bool debug) 
  : dft_BasicLinProbMgr(numUnknownsPerNode, solverOptions, solverParams, comm),
    isLinear_(false),
    debug_(debug),
    curRowA12_(-1),
    curRowA21_(-1) {

  return;
  }*/
//=============================================================================
dft_PolyLinProbMgr::dft_PolyLinProbMgr(int numUnknownsPerNode, Teuchos::ParameterList * parameterList, MPI_Comm comm, bool debug) 
  : dft_BasicLinProbMgr(numUnknownsPerNode, parameterList, comm),
    isLinear_(false),
    debug_(false),
    hasPoisson_(false),
    curRowA12_(-1),
    curRowA21_(-1) {
  parameterList_->set("P_location", 0); //change
  parameterList_->set("F_location", 0); //change
  return;
  }
//=============================================================================
dft_PolyLinProbMgr::~dft_PolyLinProbMgr() {
   return;
}
//=============================================================================
int dft_PolyLinProbMgr::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning

	TEUCHOS_TEST_FOR_EXCEPTION((numGlobalNodes_==0 ||
					   numGlobalBoxNodes_==0 ||
					   gEquations_.Length()==0 ||
					   cmsEquations_.Length()==0 ||
					   densityEquations_.Length()==0), std::logic_error, 
					   "One or more set methods not called.");
  //Not checking if poissonEquations_.Length()==0 because don't HAVE to have Poisson equations      
  //Not checking if gInvEquations_.Length()==0 because we don't have to have G inv equations
  
  int poissonInA11 = Teuchos::getParameter<int>(*parameterList_, "P_location");
  int F_location = Teuchos::getParameter<int>(*parameterList_, "F_location");

  // Fill physics ordering vector with the concatenated contents of the IDs for all physics types
  // Load Schur block mappings
  physicsOrdering_.Size(numUnknownsPerNode_);
  physicsIdToSchurBlockId_.Size(numUnknownsPerNode_);
  isCmsEquation_.Size(numUnknownsPerNode_); 
  isDensityEquation_.Size(numUnknownsPerNode_);
  isPoissonEquation_.Size(numUnknownsPerNode_);
  int * ptr = physicsOrdering_.Values();
  for (int i=0; i<gEquations_.Length(); i++) {
    *ptr++ = gEquations_[i];
    physicsIdToSchurBlockId_[gEquations_[i]] = 1;
  }
  for (int i=0; i<gInvEquations_.Length(); i++) {
    *ptr++ = gInvEquations_[i];
    physicsIdToSchurBlockId_[gInvEquations_[i]] = 1;
  }
  if (poissonInA11) {
    for (int i=0; i<poissonEquations_.Length(); i++) {
      *ptr++ = poissonEquations_[i];
      physicsIdToSchurBlockId_[poissonEquations_[i]] = 1; //so it's in A11
      isPoissonEquation_[poissonEquations_[i]] = 1;
    }
  }
  else {
    for (int i=0; i<poissonEquations_.Length(); i++) {
      *ptr++ = poissonEquations_[i];
      physicsIdToSchurBlockId_[poissonEquations_[i]] = 2; //so it's in A22
      isPoissonEquation_[poissonEquations_[i]] = 1;
    }
  }
  if (F_location == 1) { //F in NE
    for (int i=0; i<cmsEquations_.Length(); i++) {
      *ptr++ = cmsEquations_[i];
      physicsIdToSchurBlockId_[cmsEquations_[i]] = 2;
      isCmsEquation_[cmsEquations_[i]] = 1;
    }
    for (int i=0; i<densityEquations_.Length(); i++) {
      *ptr++ = densityEquations_[i];
      physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
      isDensityEquation_[densityEquations_[i]] = 1;
    }
  }
  else { //F in SW
    for (int i=0; i<densityEquations_.Length(); i++) {
      *ptr++ = densityEquations_[i];
      physicsIdToSchurBlockId_[densityEquations_[i]] = 2;
      isDensityEquation_[densityEquations_[i]] = 1;
    }
    for (int i=0;  i<cmsEquations_.Length(); i++) {
      *ptr++ = cmsEquations_[i];
      physicsIdToSchurBlockId_[cmsEquations_[i]] = 2;
      isCmsEquation_[cmsEquations_[i]] = 1;
    }
  }

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.Size(numUnknownsPerNode_);
  for (int i=0; i<physicsOrdering_.Length(); i++) solverOrdering_[physicsOrdering_[i]]=i;

  const int numUnks = numOwnedNodes_*numUnknownsPerNode_;
  const int numUnks1 = numOwnedNodes_*(gEquations_.Length()+gInvEquations_.Length());
  const int numUnks2 = numOwnedNodes_*(cmsEquations_.Length()+densityEquations_.Length());
  const int numUnksP = numOwnedNodes_*(poissonEquations_.Length());
  if (numUnksP > 0) hasPoisson_ = true;
  assert(numUnks==(numUnks1+numUnks2+numUnksP)); //Sanity test
  const int numCms = numOwnedNodes_*(cmsEquations_.Length());
  const int numDensity = numOwnedNodes_*(densityEquations_.Length());
  Epetra_IntSerialDenseVector globalGIDList(numUnks);

  int * GIDs = ownedMap_->MyGlobalElements();
  int k=0;
  for (int i=0; i<numUnknownsPerNode_; i++) {
    int ii=physicsOrdering_[i];
    for (int j=0; j<numOwnedNodes_; j++) 
      globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
  }

  ptr = globalGIDList.Values();
  globalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks, ptr, 0, comm_));
  if (poissonInA11) {
    block1RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks1+numUnksP, ptr, 0, comm_));
    block2RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks2, ptr+numUnks1+numUnksP, 0, comm_));
  }
  else {
    block1RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks1, ptr, 0, comm_));
    block2RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks2+numUnksP, ptr+numUnks1, 0, comm_));
  }
  
  if (F_location == 1) { //F in NE
    cmsRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numCms, ptr+numUnks1+numUnksP, 0, comm_));
    densityRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numDensity, ptr+numUnks1+numUnksP+numCms, 0, comm_));
  }
  else { //F in SW
    densityRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numDensity, ptr+numUnks1+numUnksP, 0, comm_));
    cmsRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numCms, ptr+numUnks1+numUnksP+numDensity, 0, comm_));
  }
  if (hasPoisson_) {
    poissonRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnksP, ptr+numUnks1, 0, comm_));
    if (poissonInA11) {
      extraRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks1, ptr, 0, comm_));
      A11_ = Teuchos::rcp(new dft_PolyA11_Coulomb_Epetra_Operator(*ownedMap_, *block1RowMap_, *extraRowMap_, *poissonRowMap_, parameterList_));
      // A11_ = Teuchos::rcp(new dft_PolyA11_Coulomb_Epetra_Operator(*ownedMap_, *block1RowMap_, *extraRowMap_, *poissonRowMap_, solverOptions_, solverParams_));
    }
    else { //poisson in A22
      extraRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks2, ptr+numUnks1+numUnksP, 0, comm_));
      A11_ = Teuchos::rcp(new dft_PolyA11_Epetra_Operator(*ownedMap_, *block1RowMap_));
    }
  }
  else { //does not have Poisson equations
    poissonRowMap_ = Teuchos::null;
    extraRowMap_ = Teuchos::null;
    A11_ = Teuchos::rcp(new dft_PolyA11_Epetra_Operator(*ownedMap_, *block1RowMap_));
  }
  
  A12_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *block1RowMap_, 0)); A12_->SetLabel("PolyLinProbMgr::A12");
  A21_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *block2RowMap_, 0)); A21_->SetLabel("PolyLinProbMgr::A21");

  if (hasPoisson_ && !poissonInA11) {
    A22_ = Teuchos::rcp(new dft_PolyA22_Coulomb_Epetra_Operator(*cmsRowMap_, *densityRowMap_, *poissonRowMap_, *extraRowMap_, *block2RowMap_, parameterList_));
    //A22_ = Teuchos::rcp(new dft_PolyA22_Coulomb_Epetra_Operator(*cmsRowMap_, *densityRowMap_, *poissonRowMap_, *extraRowMap_, *block2RowMap_, solverOptions_, solverParams_));
  } else {
    A22_ = Teuchos::rcp(new dft_PolyA22_Epetra_Operator(*cmsRowMap_, *densityRowMap_, *block2RowMap_, parameterList_));
    //A22_ = Teuchos::rcp(new dft_PolyA22_Epetra_Operator(*cmsRowMap_, *densityRowMap_, *block2RowMap_, solverOptions_, solverParams_));
  }
  if (debug_) {
    globalMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalRowMap_, 0));
    globalMatrix_->SetLabel("PolyLinProbMgr::globalMatrix");
  }
  else
    globalMatrix_ = Teuchos::null; // not used by this solver
  
  globalRhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalLhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));

  int offset2 = numUnks1;

  if (hasPoisson_ && poissonInA11) {
    offset2 += numUnksP;
  }
  rhs1_ = Teuchos::rcp(new Epetra_Vector(View, *block1RowMap_, globalRhs_->Values()));
  rhs2_ = Teuchos::rcp(new Epetra_Vector(View, *block2RowMap_, globalRhs_->Values()+offset2));
  rhsSchur_ = Teuchos::rcp(new Epetra_Vector(*rhs2_));
  lhs1_ = Teuchos::rcp(new Epetra_Vector(View, *block1RowMap_, globalLhs_->Values()));
  lhs2_ = Teuchos::rcp(new Epetra_Vector(View, *block2RowMap_, globalLhs_->Values()+offset2));

  schurOperator_ = Teuchos::rcp(new dft_Schur_Epetra_Operator(A11_.get(), A12_.get(), A21_.get(), A22_.get()));
  implicitProblem_ = Teuchos::rcp(new Epetra_LinearProblem(schurOperator_.get(), lhs2_.get(), rhsSchur_.get()));

  ownedToBoxImporter_ = Teuchos::rcp(new Epetra_Import(*boxMap_, *ownedMap_));

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::initializeProblemValues() {
  
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
  A22_->setFieldOnDensityIsLinear(isLinear_);  // Set current state of linearity for F
  A22_->initializeProblemValues();
  
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value) {

  int schurBlockRow = physicsIdToSchurBlockId_[ownedPhysicsID];
  int schurBlockCol = physicsIdToSchurBlockId_[boxPhysicsID];
  int rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  int colGID = boxToSolverGID(boxPhysicsID, boxNode);

  //cout << std::setprecision(2);
  //cout << "A[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode
  //     << "][boxPhysicsID="  <<boxPhysicsID  <<"][boxNode="  <<boxNode
  //     << "][rowGID="        <<rowGID        <<"][colGID="   <<colGID  
  //     << "] = " << value << endl;

  if (schurBlockRow==1 && schurBlockCol==1) { // A11 block
    A11_->insertMatrixValue(solverOrdering_[ownedPhysicsID], ownedMap_->GID(ownedNode), rowGID, colGID, value); 
  }
  else if (schurBlockRow==2 && schurBlockCol==2) { // A22 block
    // if poisson then blockColFlag = 0
    // if density then blockColFlag = 1
    // if cms then blockColFlag = 2
    if (isCmsEquation_[boxPhysicsID]) 
      A22_->insertMatrixValue(rowGID, colGID, value, 2); 
    else if (isDensityEquation_[boxPhysicsID])
      A22_->insertMatrixValue(rowGID, colGID, value, 1); 
    else if (isPoissonEquation_[boxPhysicsID])
      A22_->insertMatrixValue(rowGID, colGID, value, 0); 
    else 
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, "Unknown box physics ID in A22.");
  }
  else if (schurBlockRow==2 && schurBlockCol==1) { // A21 block
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
int dft_PolyLinProbMgr::insertRowA12() {
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
int dft_PolyLinProbMgr::insertRowA21() {
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
int dft_PolyLinProbMgr::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) {
    insertRowA12(); // Dump any remaining entries
    A12_->FillComplete(*block2RowMap_,*block1RowMap_);
    A12_->OptimizeStorage();
    insertRowA21(); // Dump any remaining entries
    A21_->FillComplete(*block1RowMap_,*block2RowMap_);
    A21_->OptimizeStorage();

    if (debug_) globalMatrix_->FillComplete();
    //if (debug_) globalMatrix_->OptimizeStorage(); //should this be included?

  }
  //std::cout << *A12_ << endl 
  //          << *A21_ << endl;

  A11_->finalizeProblemValues();
  A22_->finalizeProblemValues();

  //cout << "Inf Norm of A12 = " << A12_->NormInf() << endl;
  //cout << "Inf Norm of A21 = " << A21_->NormInf() << endl;
  //cout << "Fro Norm of A12 = " << A12_->NormFrobenius() << endl;
  //cout << "Fro Norm of A21 = " << A21_->NormFrobenius() << endl;
  //Check(true);

  //cout << *A21_ << endl;

  //Check(true);
  isLinearProblemSet_ = true;
  firstTime_ = false;

  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::setupSolver() {

	TEUCHOS_TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error, 
					   "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");
		
	
  schurOperator_->ComputeRHS(*rhs1_, *rhs2_, *rhsSchur_);
  if (solver_ != Teuchos::null) return(0); //All done if solver already defined.
	
  solver_ = Teuchos::rcp(new AztecOO(*implicitProblem_));
  solver_->SetParameters(*parameterList_);
  //if (solverOptions_!=0) solver_->SetAllAztecOptions(solverOptions_);
  // if (solverParams_!=0) solver_->SetAllAztecParams(solverParams_);

  //const int * options = solver_->GetAllAztecOptions();
  //const double * params = solver_->GetAllAztecParams();

  solver_->SetAztecOption(AZ_scaling, AZ_none); 
  int maxiter = 500;
  solver_->SetAztecOption(AZ_max_iter, maxiter);
  solver_->SetAztecOption(AZ_kspace, maxiter); 
  //solver_->SetAztecParam(AZ_tol, 1.0e-12); 
  //solver_->SetAztecOption(AZ_conv, AZ_noscaled); 
  solver_->SetPrecOperator(A22_.get());
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
int dft_PolyLinProbMgr::solve() {

  //writeMatrix("2D.mm", "Small Polymer Matrix", "Global Matrix from Small Polymer Problem");
  //abort();

  // const int * options = solver_->GetAllAztecOptions();
  // const double * params = solver_->GetAllAztecParams();
  // solver_->Iterate(options[AZ_max_iter], params[AZ_tol]); // Try to solve

  solver_->Iterate(Teuchos::getParameter<int>(*parameterList_, "Max_iter"), Teuchos::getParameter<double>(*parameterList_, "Tol")); // Try to solve

  schurOperator_->ComputeX1(*rhs1_, *lhs2_, *lhs1_); // Compute rest of solution

  //cout << "X1 = " << *lhs1_ << endl;
  //cout << "X2 = " << *lhs2_ << endl;

  if (debug_) {
    Epetra_Vector tmpRhs(*globalRowMap_);
    Epetra_Vector tmprhs1(View, *block1RowMap_, tmpRhs.Values());
    Epetra_Vector tmprhs2(View, *block2RowMap_, tmpRhs.Values()+block1RowMap_->NumMyElements());
    
    schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, tmprhs1, tmprhs2);
    
    tmpRhs.Update(-1.0, *globalRhs_, 1.0);
    double resid=0.0;
    tmpRhs.Norm2(&resid);
    std::cout << "Global Residual for solution = " << resid << std::endl;
    bool writeMatrixNow = false;
    if (writeMatrixNow) {
      writeMatrix("A.dat", "GlobalMatrix", "GlobalMatrix");
      writeLhs("x.dat");
      writeRhs("b.dat");
      writePermutation("p.dat");
      //abort();
    }
  }

  //solver_->AdaptiveIterate(solverOptions_[AZ_max_iter], 5, solverParams_[AZ_tol]); // Try to solve

  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::applyMatrix(const double** x, double** b) const {
  
  setLhs(x);
  schurOperator_->ApplyGlobal(*lhs1_, *lhs2_, *rhs1_, *rhs2_);
  getRhs(b);

  return(0);
}
//=============================================================================
  int dft_PolyLinProbMgr::Check(bool verbose) const {

  int ierr1 = A11_->Check(verbose);
  int ierr2 = A22_->Check(verbose);
  if (ierr1!=0 || ierr2!=0) return(-1);
  return(0);
    
  }
//=============================================================================
int dft_PolyLinProbMgr::writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const  {
cout << "in poly version of writeMatrix filename="<<filename<<endl;
/*  if (debug_)*/
    return(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_, matrixName, matrixDescription));
/*  else
    return(-1);*/ // Not available if not in debug mode
}
