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
#include "Epetra_MpiComm.h"
#include "AztecOO.h"


//=============================================================================
dft_SolverManager::dft_SolverManager(int numUnknownsPerNode, int * unknownToPhysicsType, int * solverOptions, double * solverParams, MPI_Comm comm) 
  : numUnknownsPerNode_(numUnknownsPerNode),
    numOwnedNodes_(0),
    numBoxNodes_(0),
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
    ownedMap_(0),
    boxMap_(0),
    ownedToBoxImporter_(0),
    globalMatrix_(0),
    globalRhs_(0),
    globalLhs_(0),
    solver_(0),
    isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    groupByPhysics_(true) {

  return;
}
//=============================================================================
dft_SolverManager::~dft_SolverManager() {
   return;
}
//=============================================================================
int dft_SolverManager::setNodalRowMap(int numOwnedNodes, int * GIDs, int nx=0, int ny = 1, int nz = 1) {
  numOwnedNodes_ = numOwnedNodes;

  const int numUnks = numOwnedNodes*numUnknownsPerNode_
  Epetra_IntSerialDenseVector globalGIDList(numUnks);

  int k=0;
  if (groupByPhysics_) 
    for (int i=0; i<numUnknownsPerNode_; i++)
      for (int j=0; j<numOwnedNodes_; j++) globalGIDList[k++] = i + GIDs[j]*numOwnedNodes;
  else
    for (int i=0; i<numUnknownsPerNode_; i++)
      for (int j=0; j<numOwnedNodes_; j++) globalGIDList[k++] = i*numOwnedNodes + GIDs[j];

  globalRowMap_ = Teuchos::rcp(new Epetra_Map(-1, numUnks, globalGIDList.Values(), 0, comm_));

  ownedMap_ = Teuchos::rcp(new Epetra_Map(-1, numOwnedNodes, GIDs, 0, comm_));
  return(0);
}
//=============================================================================
int dft_SolverManager::setNodalColMap(int numBoxNodes, int * GIDs, int nx=0, int ny = 1, int nz = 1) {
  
  numBoxNodes_ = numBoxNodes;

  boxMap_ = Teuchos::rcp(new Epetra_Map(-1, numBoxNodes, GIDs, 0, comm_));

  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning
  
  globalMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *globalRowMap_, 0));
  globalRhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalLhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalProblem_ = Teuchos::rcp(new Epetra_LinearProblem(globalMatrix_, globalLhs_, globalRhs_));
    
  ownedToBoxImporter_ = Teuchos::rcp(new Epetra_Import(*boxMap_, *ownedMap_));

  isBlockStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_SolverManager::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime) {
    globalMatrix_->PutScalar(0.0);
    globalRhs_->PutScalar(0.0);
    globalLhs_->PutScalara(0.0);
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
  if (firstTime)
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

  isLinearProblemSet_ = true;
  return(0);
}
//=============================================================================
int dft_SolverManager::setBlockMatrixReadOnly(int rowPhysicsID, int colPhysicsID, bool readOnly) {
  
  return(-1);  // Not implemented yet.
}
//=============================================================================
int dft_SolverManager::setRhs(const double ** b) {

  for (int i=0; i<numUnknownsPerNode_; i++)
    for (int j=0; j<numOwnedNodes_; j++)
      globalRhs[ownedToSolverLID(i,j)] = b[i][j];
  
  return(0);
}
//=============================================================================
int dft_SolverManager::setLhs(double ** x) const {

  Epetra_SerialDenseVector xtmp(numOwnedNodes_); // Temp vector to hold local x values
  for (int i=0; i<numUnknownsPerNode_; i++) {
    exportC2R(x[i], xtmp.Values()); // Use simple import
    for (int j=0; j<numOwnedNodes_; j++)
      globalLhs[ownedToSolverLID(i,j)] = xtmp[j];
  }
  return(0);
}
//=============================================================================
int dft_SolverManager::getLhs(double ** x) const {

  Epetra_SerialDenseVector xtmp(numOwnedNodes_); // Temp vector to hold local x values
  for (int i=0; i<numUnknownsPerNode_; i++) {
    for (int j=0; j<numOwnedNodes_; j++)
      xtmp[j] = globalLhs[ownedToSolverLID(i,j)];
    importR2C(xtmp.Values(), x[i]); // Use simple import
  }
  return(0);
}
//=============================================================================
int dft_SolverManager::getRhs(double ** b) const {

  for (int i=0; i<numUnknownsPerNode_; i++)
    for (int j=0; j<numOwnedNodes_; j++)
      b[i][j] = globalRhs[ownedToSolverLID(i,j)];
  
  
  return(0);
}
//=============================================================================
int dft_SolverManager::setupSolver() {

  solver_ = Teuchos::rcp(new AztecOO(globalProblem_));
  solver_->SetAztecOption(AZ_solver, AZ_gmres);
  solver_->SetAztecOption(AZ_precond, AZ_dom_decomp);
  solver_->SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  solver_->SetAztecParam(AZ_ilut_fill, 4.0);
  solver_->SetAztecParam(AZ_drop, 0.0);

  return(0);
}
//=============================================================================
int dft_SolverManager::solve() {
  
  solver_->Iterate(200, 1.0E-8); // Try to solve
  return(0);
}
//=============================================================================
int dft_SolverManager::applyMatrix(const double** x, double** b) const {
  
  setLhs(x);
  globalMatrix_->Apply(false, *globalLhs_, *globalRhs_);
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
  
  Epetra_Vector owned(View, ownedMap, (double *) aOwned);
  Epetra_Vector box(View, boxMap, aBox);
  
  box.Import(owned, *ownedToBoxImporter_, Insert);

  return(0);
}
//=============================================================================
int dft_SolverManager::exportC2R(const double* aBox, double* aOwned) const {
  
  Epetra_Vector owned(View, ownedMap, aOwned);
  Epetra_Vector box(View, boxMap, (double *) aBox);
  
  owned.Export(box, *ownedToBoxImporter_, Zero); // Use importer, but zero out off-processor contributions.

  return(0);
}
//=============================================================================

