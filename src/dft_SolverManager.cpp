/*@HEADER
// ***********************************************************************
// 
//                Tramonto: Molecular Theories Modeling Code
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// Questions? Contact Laura J.D. Frink (ljfrink@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/

#include "dft_SolverManager.hpp"
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


//=============================================================================
dft_SolverManager::dft_SolverManager(int numUnknownsPerNode, int * unknownToPhysicsType, int * solverOptions, double * solverParams, MPI_Comm comm) 
  : numUnknownsPerNode_(numUnknownsPerNode),
    solverOptions_(solverOptions),
    solverParams_(solverParams),
    numOwnedNodes_(0),
    numBoxNodes_(0),
    numGlobalNodes_(0),
    numMatrixBlocks_(0),
    comm_(Epetra_MpiComm(comm)),
    blockMatrix_(0),
    blockMatrixReadOnly_(0),
    blockMatrixIsVbr_(0),
    blockGraph_(0),
    blockLhs_(0),
    blockRhs_(0),
    rowMaps_(0),
    colMaps_(0),
    isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    groupByPhysics_(true),
    firstTime_(true) {

  return;
}
//=============================================================================
dft_SolverManager::~dft_SolverManager() {
   return;
}
//=============================================================================
int dft_SolverManager::setNodalRowMap(int numOwnedNodes, int * GIDs, int nx, int ny, int nz) {
  numOwnedNodes_ = numOwnedNodes;
  comm_.SumAll(&numOwnedNodes_, &numGlobalNodes_, 1);

  ownedMap_ = Teuchos::rcp(new Epetra_Map(-1, numOwnedNodes, GIDs, 0, comm_));
  //std::cout << " Owned Map" << *ownedMap_.get() << std::endl;
  return(0);
}
//=============================================================================
int dft_SolverManager::setNodalColMap(int numBoxNodes, int * GIDs, int nx, int ny, int nz) {
  
  numBoxNodes_ = numBoxNodes;

  boxMap_ = Teuchos::rcp(new Epetra_Map(-1, numBoxNodes, GIDs, 0, comm_));
  //std::cout << " Box Map" << *boxMap_.get() << std::endl;

  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning
  
  const int numUnks = numOwnedNodes_*numUnknownsPerNode_;
  Epetra_IntSerialDenseVector globalGIDList(numUnks);

  int * GIDs = ownedMap_->MyGlobalElements();
  int k=0;
  if (groupByPhysics_) 
    for (int i=0; i<numUnknownsPerNode_; i++)
      for (int j=0; j<numOwnedNodes_; j++) 
	globalGIDList[k++] = i*numGlobalNodes_ + GIDs[j];
  else
    for (int j=0; j<numOwnedNodes_; j++) 
      for (int i=0; i<numUnknownsPerNode_; i++)
	globalGIDList[k++] = i + GIDs[j]*numUnknownsPerNode_;

  globalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks, globalGIDList.Values(), 0, comm_));

  //std::cout << " Global Row Map" << *globalRowMap_.get() << std::endl;

  globalMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalRowMap_, 0));
  globalRhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalLhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalProblem_ = Teuchos::rcp(new Epetra_LinearProblem(globalMatrix_.get(), globalLhs_.get(), globalRhs_.get()));
    
  ownedToBoxImporter_ = Teuchos::rcp(new Epetra_Import(*(boxMap_.get()), *(ownedMap_.get())));

  isBlockStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_SolverManager::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    globalMatrix_->PutScalar(0.0);
    globalRhs_->PutScalar(0.0);
    globalLhs_->PutScalar(0.0);
  }
  
  return(0);
}
//=============================================================================
int dft_SolverManager::insertRhsValue(int ownedPhysicsID, int ownedNode, double value) {

  int rhsLID = ownedToSolverLID(ownedPhysicsID, ownedNode); // Get solver LID
  (*globalRhs_)[rhsLID] = value;
  return(0);
}
//=============================================================================
int dft_SolverManager::insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value) {

  int rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  int colGID = boxToSolverGID(boxPhysicsID, boxNode);
  if (firstTime_)
    globalMatrix_->InsertGlobalValues(rowGID, 1, &value, &colGID);
  else
    globalMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
  
  return(0);
}
//=============================================================================
int dft_SolverManager::insertMatrixValues(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int * boxNodeList, double * values, int numEntries) {
  
  for (int i=0; i<numEntries; i++) insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNodeList[i], values[i]);

  return(0);
}
//=============================================================================
int dft_SolverManager::insertMatrixValues(int ownedPhysicsID, int ownedNode, int * boxPhysicsIDList, int boxNode, double * values, int numEntries) {
  
  for (int i=0; i<numEntries; i++) insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsIDList[i], boxNode, values[i]);

  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  globalMatrix_->FillComplete();
  globalMatrix_->OptimizeStorage();

  //std::cout << *globalMatrix_.get();

  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//=============================================================================
int dft_SolverManager::setBlockMatrixReadOnly(int rowPhysicsID, int colPhysicsID, bool readOnly) {
  
  return(-1);  // Not implemented yet.
}
//=============================================================================
int dft_SolverManager::setRhs(const double ** b) {

  double * tmp = globalRhs_->Values();
  for (int i=0; i<numUnknownsPerNode_; i++)
    for (int j=0; j<numOwnedNodes_; j++)
      tmp[ownedToSolverLID(i,j)] = b[i][j];
  
  return(0);
}
//=============================================================================
int dft_SolverManager::setLhs(const double ** x) const {

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
int dft_SolverManager::getLhs(double ** x) const {

  double * tmp = globalLhs_->Values();
  Epetra_SerialDenseVector xtmp(numOwnedNodes_); // Temp vector to hold local x values
  for (int i=0; i<numUnknownsPerNode_; i++) {
    for (int j=0; j<numOwnedNodes_; j++)
      xtmp[j] = tmp[ownedToSolverLID(i,j)];
    importR2C(xtmp.Values(), x[i]); // Use simple import
  }
  return(0);
}
//=============================================================================
int dft_SolverManager::getRhs(double ** b) const {

  double * tmp = globalRhs_->Values();
  for (int i=0; i<numUnknownsPerNode_; i++)
    for (int j=0; j<numOwnedNodes_; j++)
      b[i][j] = tmp[ownedToSolverLID(i,j)];
  
  
  return(0);
}
//=============================================================================
int dft_SolverManager::setupSolver() {

  solver_ = Teuchos::rcp(new AztecOO(*(globalProblem_.get())));
  solver_->SetAllAztecOptions(solverOptions_);
  solver_->SetAllAztecParams(solverParams_);
  //solver_->SetAztecOption(AZ_solver, AZ_gmres);
  //solver_->SetAztecOption(AZ_precond, AZ_dom_decomp);
  //solver_->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  //solver_->SetAztecParam(AZ_ilut_fill, 4.0);
  //solver_->SetAztecParam(AZ_drop, 0.0);

  return(0);
}
//=============================================================================
int dft_SolverManager::solve() {
  
  solver_->Iterate(solverOptions_[AZ_max_iter], solverParams_[AZ_tol]); // Try to solve
  return(0);
}
//=============================================================================
int dft_SolverManager::applyMatrix(const double** x, double** b) const {
  
  setLhs(x);
  globalMatrix_->Apply(*globalLhs_.get(), *globalRhs_.get());
  getRhs(b);
  
  return(0);
}
//=============================================================================
int dft_SolverManager::importR2C(const double** xOwned, double** xBox) const {
  
  for (int i=0; i<numUnknownsPerNode_; i++) importR2C(xOwned[i], xBox[i]);

  return(0);
}
//=============================================================================
int dft_SolverManager::importR2C(const double* aOwned, double* aBox) const {
  
  Epetra_Vector owned(View, *ownedMap_.get(), (double *) aOwned);
  Epetra_Vector box(View, *boxMap_.get(), aBox);
  
  box.Import(owned, *ownedToBoxImporter_.get(), Insert);

  return(0);
}
//=============================================================================
int dft_SolverManager::exportC2R(const double* aBox, double* aOwned) const {
  
  Epetra_Vector owned(View, *ownedMap_.get(), aOwned);
  Epetra_Vector box(View, *boxMap_.get(), (double *) aBox);
  
  owned.Export(box, *ownedToBoxImporter_.get(), Zero); // Use importer, but zero out off-processor contributions.

  return(0);
}
//=============================================================================

