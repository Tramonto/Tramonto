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

#include "dft_PolyLinProbMgr.hpp"
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


//=============================================================================
dft_PolyLinProbMgr::dft_PolyLinProbMgr(int numUnknownsPerNode, int * solverOptions, double * solverParams, MPI_Comm comm) 
  : dft_BasicLinProbMgr(numUnknownsPerNode, solverOptions, solverParams, comm) {
  

  return;
}
//=============================================================================
dft_PolyLinProbMgr::~dft_PolyLinProbMgr() {
   return;
}
//=============================================================================
int dft_PolyLinProbMgr::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning

  if (numGlobalNodes_==0 ||
      numGlobalBoxNodes_==0 ||
      gInvEquations_.Length()==0 ||
      gInvEquations_.Length()==0 ||
      cmsEquations_.Length()==0 ||
      densityEquations_.Length()==0) return(-1); // Error: One or more set methods not called
      
  
  // Fill physics ordering vector with the concatenated contents of the IDs for all physics types

  physicsOrdering_.Size(numUnknownsPerNode_);
  int * ptr = physicsOrdering_.Values();
  for (int i=0; i<gEquations_.Length(); i++) *ptr++ = gEquations_[i];
  for (int i=0; i<gInvEquations_.Length(); i++) *ptr++ = gInvEquations_[i];
  for (int i=0; i<cmsEquations_.Length(); i++) *ptr++ = cmsEquations_[i];
  for (int i=0; i<densityEquations_.Length(); i++) *ptr++ = densityEquations_[i];

  const int numUnks = numOwnedNodes_*numUnknownsPerNode_;
  const int numUnks1 = gEquations_.Length()+gInvEquations_.Length();
  const int numUnks2 = cmsEquations_.Length()+densityInvEquations_.Length();
  assert(numUnks==(numUnks1+numUnks2));  // Sanity test
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
  block1RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numOwnedNodes_*numUnks1, ptr, 0, comm_));
  block2RowMap_ = Teuchos::rcp(new Epetra_Map(-1, numOwnedNodes_*numUnks2, ptr+numUnks1, 0, comm_));

  //std::cout << " Global Row Map" << *globalRowMap_.get() << std::endl;

  A11_ = Teuchos::rcp(new dft_PolyA11_Epetra_Operator(ownedMap_, gEquations.Length()));
  A12_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *block1RowMap_, 0));
  A21_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *block2RowMap_, 0));
  A22_ = Teuchos::rcp(new dft_PolyA22_Epetra_Operator(ownedMap_, cmsEquations_.Length()+densityEquations_.Length()));
  globalMatrix = null; // not used by this solver

  globalRhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalLhs_ = Teuchos::rcp(new Epetra_Vector(*globalRowMap_));
  globalProblem_ = Teuchos::rcp(new Epetra_LinearProblem(globalMatrix_.get(), globalLhs_.get(), globalRhs_.get()));
    
  ownedToBoxImporter_ = Teuchos::rcp(new Epetra_Import(*(boxMap_.get()), *(ownedMap_.get())));

  isBlockStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::initializeProblemValues() {
  
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
int dft_PolyLinProbMgr::insertRhsValue(int ownedPhysicsID, int ownedNode, double value) {

  int rhsLID = ownedToSolverLID(ownedPhysicsID, ownedNode); // Get solver LID
  (*globalRhs_)[rhsLID] += value;
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value) {

  int rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  int colGID = boxToSolverGID(boxPhysicsID, boxNode);
  if (firstTime_)
    globalMatrix_->InsertGlobalValues(rowGID, 1, &value, &colGID);
  else
    globalMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
  
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  globalMatrix_->FillComplete();
  globalMatrix_->OptimizeStorage();

  //std::cout << *globalMatrix_.get();

  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::setupSolver() {

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
int dft_PolyLinProbMgr::solve() {
  
  //writeMatrix("2D.mm", "Small Polymer Matrix", "Global Matrix from Small Polymer Problem");
  //abort();
  solver_->Iterate(solverOptions_[AZ_max_iter], solverParams_[AZ_tol]); // Try to solve
  //solver_->AdaptiveIterate(solverOptions_[AZ_max_iter], 5, solverParams_[AZ_tol]); // Try to solve
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::applyMatrix(const double** x, double** b) const {
  
  setLhs(x);
  globalMatrix_->Apply(*globalLhs_.get(), *globalRhs_.get());
  getRhs(b);
  
  return(0);
}
//=============================================================================
int dft_PolyLinProbMgr::writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const  {
    return(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_.get(), matrixName, matrixDescription));
}
//=============================================================================

