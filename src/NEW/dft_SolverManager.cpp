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
#include "Epetra_MpiComm.h"


//=============================================================================
dft_SolverManager::dft_SolverManager(int numUnknownsPerNode, int * unknownToPhysicsType, int * solverOptions, double * solverParams, MPI_Comm comm) 
  : numUnknownsPerNode_(numUnknownsPerNode),
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
    isLinearProblemSet_(false) {

  return;
}
//=============================================================================
dft_SolverManager::~dft_SolverManager() {
 
  return;
}
//=============================================================================
int dft_SolverManager::setNodalRowMap(int numOwnedNodes, int * GIDs, int nx=0, int ny = 1, int nz = 1) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::setNodalColMap(int numBoxNodes, int * GIDs, int nx=0, int ny = 1, int nz = 1) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeBlockStructure() {

  if (isBlockStructureSet_) return(1); // Already been here, return warning
  
  bool useDefaultSolver = true;  // Make this always true for now.
  
  // Compute owned Physics/Node to solver GID translation mappings, then the same for box physics/nodes.
  compute
  if (useDefaultSolver) {
    isBlockStructureSet_ = true;
  return(0);
}
//=============================================================================
int dft_SolverManager::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem
  
  return(0);
}
//=============================================================================
int dft_SolverManager::insertRhsValue(int ownedPhysicsID, int ownedNode, double value) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::insertMatrixValues(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int * boxNodeList, double * values, int numEntries) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::insertMatrixValues(int ownedPhysicsID, int ownedNode, int * boxPhysicsIDList, int boxNode, double * values, int numEntries) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do
  isLinearProblemSet_ = true;
  return(0);
}
//=============================================================================
int dft_SolverManager::setBlockMatrixReadOnly(int rowPhysicsID, int colPhysicsID, bool readOnly) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::setRhs(const double ** b) {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::getLhs(double ** x) const {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::getRhs(double ** b) const {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::setupSolver() {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::solve() {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::applyMatrix(const double** x, double** b) const {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::importR2C(const double** xOwned, double** xBox) const {
  
  return(0);
}
//=============================================================================
int dft_SolverManager::importR2C(const double* aOwned, double* aBox) const {
  
  return(0);
}
//=============================================================================

