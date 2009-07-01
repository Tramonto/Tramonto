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

#include "dft_BasicLinProbMgr.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_BlockMap.h"
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


//=============================================================================
/*dft_BasicLinProbMgr::dft_BasicLinProbMgr(int numUnknownsPerNode, int * solverOptions, double * solverParams, MPI_Comm comm) 
  : numUnknownsPerNode_(numUnknownsPerNode),
    solverOptions_(solverOptions),
    solverParams_(solverParams),
    numOwnedNodes_(0),
    numBoxNodes_(0),
    numGlobalNodes_(0),
    numGlobalBoxNodes_(0),
    numCoarsenedNodes_(0),
    numGlobalCoarsenedNodes_(0),
    comm_(Epetra_MpiComm(comm)),
    isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    groupByPhysics_(true),
    firstTime_(true),
    curRow_(-1) {

  //int tmp;
  //if (comm_.MyPID()==0) std::cin>>tmp;
  //comm_.Barrier();
  return;
  }*/
//=============================================================================
dft_BasicLinProbMgr::dft_BasicLinProbMgr(int numUnknownsPerNode, Teuchos::ParameterList * parameterList, MPI_Comm comm)
  : numUnknownsPerNode_(numUnknownsPerNode),
    parameterList_(parameterList),
    numOwnedNodes_(0),
    numBoxNodes_(0),
    numGlobalNodes_(0),
    numGlobalBoxNodes_(0),
    numCoarsenedNodes_(0),
    numGlobalCoarsenedNodes_(0),
    comm_(Epetra_MpiComm(comm)),
    isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    groupByPhysics_(true),
    firstTime_(true),
    curRow_(-1) { 
    return;
    }
//=============================================================================
dft_BasicLinProbMgr::~dft_BasicLinProbMgr() {
   return;
}
//=============================================================================
int dft_BasicLinProbMgr::setNodalRowMap(int numOwnedNodes, int * GIDs, int nx, int ny, int nz) {
  if (numGlobalNodes_!=0) return(0); // Already been here
  numOwnedNodes_ = numOwnedNodes;
  comm_.SumAll(&numOwnedNodes_, &numGlobalNodes_, 1);

  ownedMap_ = Teuchos::rcp(new Epetra_Map(-1, numOwnedNodes, GIDs, 0, comm_));
  //std::cout << " Owned Map" << *ownedMap_ << std::endl;
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::setNodalColMap(int numBoxNodes, int * GIDs, int nx, int ny, int nz) {
  
  if (numGlobalBoxNodes_!=0) return(0); // Already been here

  numBoxNodes_ = numBoxNodes;
  comm_.SumAll(&numBoxNodes_, &numGlobalBoxNodes_, 1);

  boxMap_ = Teuchos::rcp(new Epetra_Map(-1, numBoxNodes, GIDs, 0, comm_));
  //std::cout << " Box Map" << *boxMap_ << std::endl;

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::setCoarsenedNodesList(int numCoarsenedNodes, int * GIDs) {
  
  if (numGlobalCoarsenedNodes_!=0) return(0); // Already been here

  numCoarsenedNodes_ = numCoarsenedNodes;
  comm_.SumAll(&numCoarsenedNodes_, &numGlobalCoarsenedNodes_, 1);

  coarsenedNodesMap_ = Teuchos::rcp(new Epetra_Map(-1, numCoarsenedNodes, GIDs, 0, comm_));
  //std::cout << " Coarsened Nodes Map" << *coarsenedNodesMap_ << std::endl;

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning
  
  const int numUnks = numOwnedNodes_*numUnknownsPerNode_;
  Epetra_IntSerialDenseVector globalGIDList(numUnks);

  // Physics ordering for Basic Linear Problem is natural ordering:
  physicsOrdering_.Size(numUnknownsPerNode_);
  int * ptr = physicsOrdering_.Values();
  for (int i=0; i<physicsOrdering_.Length(); i++) *ptr++ = i;

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.Size(numUnknownsPerNode_);
  for (int i=0; i<physicsOrdering_.Length(); i++) solverOrdering_[physicsOrdering_[i]]=i;

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  int * GIDs = ownedMap_->MyGlobalElements();
  int k=0;
  if (groupByPhysics_) 
    for (int i=0; i<numUnknownsPerNode_; i++) {
	int ii=physicsOrdering_[i];
	for (int j=0; j<numOwnedNodes_; j++) 
	  globalGIDList[k++] = ii*numGlobalNodes_ + GIDs[j];
    }
  else
    for (int j=0; j<numOwnedNodes_; j++)
      for (int i=0; i<numUnknownsPerNode_; i++) {
	int ii=physicsOrdering_[i];
	globalGIDList[k++] = ii + GIDs[j]*numUnknownsPerNode_;
      }
  
  globalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks, globalGIDList.Values(), 0, comm_));

  //std::cout << " Global Row Map" << *globalRowMap_ << std::endl;

  globalMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalRowMap_, 0));
  globalMatrix_->SetLabel("BasicLinProbMgr::globalMatrix");
  globalRhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalLhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  implicitProblem_ = Teuchos::rcp(new Epetra_LinearProblem(globalMatrix_.get(), globalLhs_.get(), globalRhs_.get()));
    
  ownedToBoxImporter_ = Teuchos::rcp(new Epetra_Import(*boxMap_, *ownedMap_));

  isBlockStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

 // AGS: I found that we needed to initialize the matrix even the 
 // first time a matrix was filled, because matrix entries are being put
 // in on residual-only fills, which can occur before matrix fills.

 // if (!firstTime_) {
    globalMatrix_->PutScalar(0.0);
    globalRhs_->PutScalar(0.0);
    globalLhs_->PutScalar(0.0);
//  }
  
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::insertRhsValue(int ownedPhysicsID, int ownedNode, double value) {

  int rhsLID = ownedToSolverLID(ownedPhysicsID, ownedNode); // Get solver LID
  (*globalRhs_)[rhsLID] += value;
  //cout << std::setprecision(2);
  //cout << "b[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode<<"] = " << value << endl;

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value) {

  int rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  int colGID = boxToSolverGID(boxPhysicsID, boxNode);

  //cout << std::setprecision(2);
  //cout << "A[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode
  //     << "][boxPhysicsID="  <<boxPhysicsID  <<"][boxNode="  <<boxNode
  //     << "][rowGID="        <<rowGID        <<"][colGID="   <<colGID  
  //     << "] = " << value << endl;
  if (firstTime_) {
    if (rowGID!=curRow_) { 
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    }
    curRowValues_[colGID] += value;
  }
  else
    globalMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
  
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::insertRow() {
  if (curRowValues_.empty()) return(0);
  int numEntries = curRowValues_.size();
  if (numEntries>indices_.Length()) {
    indices_.Resize(numEntries);
    values_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  globalMatrix_->InsertGlobalValues(curRow_, numEntries, values_.Values(), indices_.Values());

  curRowValues_.clear();
  return(0);
}
//=============================================================================
double dft_BasicLinProbMgr::getMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode) {

  if (globalMatrix_.get()==0) {
    std::cerr << "Global Matrix is not constructed, must set debug flag to enable this feature." << std::endl;
    abort();
  }
  
  int rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  int colGID = boxToSolverGID(boxPhysicsID, boxNode);
  int numEntries;
  int * indices;
  double * values;
  if (globalMatrix_->IndicesAreGlobal()) {
    globalMatrix_->ExtractGlobalRowView(rowGID, numEntries, values, indices); // get view of current row
    for (int i=0; i<numEntries; i++) if (colGID==indices[i]) return(values[i]);
  }
  else {
    rowGID = globalMatrix_->LRID(rowGID); // get local row ID
    colGID = globalMatrix_->LCID(colGID); // get local column ID
    globalMatrix_->ExtractMyRowView(rowGID, numEntries, values, indices); // get view of current row
    for (int i=0; i<numEntries; i++) if (colGID==indices[i]) return(values[i]);
  }
  
  return(0.0);
}
//=============================================================================
int dft_BasicLinProbMgr::insertMatrixValues(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int * boxNodeList, double * values, int numEntries) {
  
  for (int i=0; i<numEntries; i++) insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNodeList[i], values[i]);

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::insertMatrixValues(int ownedPhysicsID, int ownedNode, int * boxPhysicsIDList, int boxNode, double * values, int numEntries) {
  
  for (int i=0; i<numEntries; i++) insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsIDList[i], boxNode, values[i]);

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  if (firstTime_) {
    insertRow();
    globalMatrix_->FillComplete();
    globalMatrix_->OptimizeStorage();
  }

  isLinearProblemSet_ = true;
  firstTime_ = false;

  //writeMatrix("basica.dat", "", "");

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::setRhs(const double ** b) {

  double * tmp = globalRhs_->Values();
  for (int i=0; i<numUnknownsPerNode_; i++)
    for (int j=0; j<numOwnedNodes_; j++)
      tmp[ownedToSolverLID(i,j)] = b[i][j];
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::setLhs(const double ** x) const {

  double * tmp = globalLhs_->Values();
  Epetra_SerialDenseVector xtmp(numOwnedNodes_); // Temp vector to hold local x values
  for (int i=0; i<numUnknownsPerNode_; i++) {
    exportC2R(x[i], xtmp.Values()); // Use simple import
    for (int j=0; j<numOwnedNodes_; j++)
      tmp[ownedToSolverLID(i,j)] = xtmp[j];
  }
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::getLhs(double ** x) const {

  double * tmp = globalLhs_->Values();
  Epetra_SerialDenseVector xtmp(numOwnedNodes_); // Temp vector to hold local x values
  for (int i=0; i<numUnknownsPerNode_; i++) {
    for (int j=0; j<numOwnedNodes_; j++)
      xtmp[j] = tmp[ownedToSolverLID(i,j)];
    importR2C(xtmp.Values(), x[i]); // Use simple import
    //cout << std::setprecision(2);
    //for (int j=0; j<numOwnedNodes_; j++) cout << "x["<<i<<"]["<<j<<"] = " << x[i][j] << endl;
  }
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::getRhs(double ** b) const {

  double * tmp = globalRhs_->Values();
  for (int i=0; i<numUnknownsPerNode_; i++)
    for (int j=0; j<numOwnedNodes_; j++)
      b[i][j] = tmp[ownedToSolverLID(i,j)];
  
  
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::setupSolver() {

  if (solver_ != Teuchos::null || directSolver_  != Teuchos::null) 
     return (0);  // Already here (AGS attempt to fix memory prob June 09);

  int solverInt = Teuchos::getParameter<int>(*parameterList_, "Solver");
  // int solverInt = solverOptions_[AZ_solver];
  char * solverName;

  switch (solverInt) {
  case AM_lapack: solverName = "Amesos_Lapack"; break;
  case AM_klu: solverName = "Amesos_Klu"; break;
  case AM_mumps: solverName = "Amesos_Mumps"; break;
  case AM_umfpack: solverName = "Amesos_Umfpack"; break;
  case AM_superlu: solverName = "Amesos_Superlu"; break;
  case AM_superludist: solverName = "Amesos_Superludist"; break;
  case AM_pardiso: solverName = "Amesos_Pardiso"; break;
  case AM_taucs: solverName = "Amesos_Taucs"; break;
  default: break;
  }

  //does this need this??
  /*  if (solverInt == AM_superludist) {  
    Teuchos::ParameterList& sublist = parameterList_->sublist("Superludist");
    sublist.set("Fact", "SamePattern_SameRowPerm");
    sublist.set("Equil", true);
    sublist.set("ColPerm", "MMD_AT_PLUS_A");
    sublist.set("RowPerm", "LargeDiag");
    sublist.set("ReplaceTinyPivot", true);
    sublist.set("IterRefine", "DOUBLE");
    }*/

  if (solverInt >= AM_lapack && solverInt <= AM_taucs) {
    Amesos Amesos_Factory;
    directSolver_ = Teuchos::rcp(Amesos_Factory.Create(solverName, *implicitProblem_)); 
    directSolver_->SetParameters(*parameterList_);
    directSolver_->SymbolicFactorization();
  } else {
    solver_ = Teuchos::rcp(new AztecOO(*implicitProblem_));
    //   if (solverOptions_!=0) solver_->SetAllAztecOptions(solverOptions_);
    //    if (solverParams_!=0) solver_->SetAllAztecParams(solverParams_);
    solver_->SetParameters(*parameterList_);
  }
  
  //const int * options = solver_->GetAllAztecOptions();
  //const double * params = solver_->GetAllAztecParams();
  //solver_->SetAztecOption(AZ_scaling, AZ_none); 
  //int maxiter = 500;
  //solver_->SetAztecOption(AZ_max_iter, maxiter);
  //solver_->SetAztecOption(AZ_kspace, maxiter); 
  //solver_->SetAztecOption(AZ_conv, AZ_noscaled); 
  //solver_->SetAztecOption(AZ_precond, AZ_none);
  
  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::solve() {

  int solverInt = Teuchos::getParameter<int>(*parameterList_, "Solver");
  //  int solverInt = solverOptions_[AZ_solver];
  if (solverInt >= AM_lapack && solverInt <= AM_taucs) {
    directSolver_->NumericFactorization();
    directSolver_->Solve();
  } else {
    //writeMatrix("POLY1_WJDC.mm", "Small Polymer Matrix with WJDC", "Global Matrix from Small Polymer Problem with WJDC");
    //writePermutation("POLY1_WJDC_Perm.mm");
    //abort();

    //   const int * options = solver_->GetAllAztecOptions();
    //   const double * params = solver_->GetAllAztecParams();
    //  solver_->Iterate(options[AZ_max_iter], params[AZ_tol]); // Try to solve
    solver_->Iterate(Teuchos::getParameter<int>(*parameterList_, "Max_iter"), Teuchos::getParameter<double>(*parameterList_, "Tol")); // Try to solve
  }

    bool writeMatrixNow = false;
    if (writeMatrixNow) {
      writeMatrix("A_ref.dat", "GlobalMatrix", "GlobalMatrix");
      writeLhs("x_ref.dat");
      writeRhs("b_ref.dat");
      writePermutation("p_ref.dat");
      //abort();
    }
  //solver_->AdaptiveIterate(solverOptions_[AZ_max_iter], 5, solverParams_[AZ_tol]); // Try to solve

    //abort();

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::applyMatrix(const double** x, double** b) const {
  
  setLhs(x);
  globalMatrix_->Apply(*globalLhs_, *globalRhs_);
  getRhs(b);

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::importR2C(const double** xOwned, double** xBox) const {
  for (int i=0; i<numUnknownsPerNode_; i++) importR2C(xOwned[i], xBox[i]);

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::importR2C(const double* aOwned, double* aBox) const {
  
  Epetra_Vector owned(View, *ownedMap_, (double *) aOwned);
  Epetra_Vector box(View, *boxMap_, aBox);
  
  box.Import(owned, *ownedToBoxImporter_, Insert);

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::exportC2R(const double* aBox, double* aOwned) const {
  
  Epetra_Vector owned(View, *ownedMap_, aOwned);
  Epetra_Vector box(View, *boxMap_, (double *) aBox);
  
  owned.Export(box, *ownedToBoxImporter_, Zero); // Use importer, but zero out off-processor contributions.

  return(0);
}
//=============================================================================
int dft_BasicLinProbMgr::writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const  {
    return(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_, matrixName, matrixDescription));
}
//=============================================================================
int dft_BasicLinProbMgr::writeLhs(const char * filename) const  {
    return(EpetraExt::MultiVectorToMatlabFile(filename, *globalLhs_));
}
//=============================================================================
int dft_BasicLinProbMgr::writeRhs(const char * filename) const  {
    return(EpetraExt::MultiVectorToMatlabFile(filename, *globalRhs_));
}
//=============================================================================
int dft_BasicLinProbMgr::writePermutation(const char * filename) const  {
  return(EpetraExt::BlockMapToMatrixMarketFile(filename, *globalRowMap_, " ", " ", false));
}
//=============================================================================
int dft_BasicLinProbMgr::checkPhysicsOrdering() const  {
  if (physicsOrdering_.Length()==0) {
    std::cerr << "No unknowns are registered with this problem manager." << std::endl;
    return (-1);
  }

  bool fatalError = false;

  int numUnks = physicsOrdering_.Length();
  Epetra_SerialDenseVector tmp(numUnks);
  for (int i=0; i<numUnks; i++) {
    int curID = physicsOrdering_[i];
    if (curID <0) {
      fatalError = true;
      std::cerr << "Invalid unknown number " << curID << " is less than 0." << std::endl;
    }
    else if (curID>=numUnks) {
      fatalError = true;
      std::cerr << "Invalid unknown number " << curID << " is greater than or equal to the number of unknowns (" << numUnks << ")." << std::endl;
    }
    else tmp(curID)++; // Increment counter for this ID (at the end each ID should appear exactly one time).
  }

  for (int i=0; i<numUnks; i++) {
    if (tmp[i]==0) {
      fatalError = true;
      std::cerr << "Unknown number " << i << " is not present and should be." << std::endl;
    }
    if (tmp[i]>1) {
      fatalError = true;
      std::cerr << "Unknown number " << i << " is present " << tmp[i] << " times and should be present only once." << std::endl;
    }
  }

  if (fatalError) throw tmp.ReportError("Fatal Error in specification of unknowns", -1);
  return(0);
}
