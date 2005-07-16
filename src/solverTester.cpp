//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Time.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Flops.h"
#include "Epetra_MpiComm.h"

#include "dft_PolyLinProbMgr.hpp"
#include <mpi.h>

//int gid(int physicsID, nodeID, numPhysics, numNodes) {return(
int main(int argc, char *argv[])
{


  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm( MPI_COMM_WORLD );

  if (comm.NumProc()>1) { cout << "Can only use a single processor for now!!" << endl; exit(1);}

  int numUnknownsPerNode = 6;
  int numOwnedNodes = 4;

  int densityequ[] = {0}; int numDensity = 1;
  int cmsequ[] = {1}; int numCms = 1;
  int gequ[] = {2, 4}; int numG = 2;
  int ginvequ[] = {5, 3}; int numGinv = 2;

  Epetra_Map nodalRowMap(-1, numOwnedNodes, 0, comm);
  int numGlobalNodes = nodalRowMap.NumGlobalElements();
  int * nodalElements = nodalRowMap.MyGlobalElements();
  int numMyUnknowns = numOwnedNodes*numUnknownsPerNode;
  Epetra_IntSerialDenseVector globalGids(numMyUnknowns);
  int * ptr = globalGids.Values();
  for (int i=0; i<numUnknownsPerNode; i++) 
    for (int j=0; j<numOwnedNodes; j++) *ptr++ = numGlobalNodes*i + j;

  Epetra_Map globalRowMap(-1, numMyUnknowns, globalGids, 0, comm);
  Epetra_CrsMatrix(Copy, globalRowMap);
  

  // For now we are on a single processor
  int numBoxNodes = numOwnedNodes;
  Epetra_Map nodalColMap(nodalRowMap);

  dft_PolyLinProbMgr mgr(numUnknownsPerNode, 0, 0, MPI_COMM_WORLD);
  mgr.setNodalRowMap(numOwnedNodes, nodalRowMap.MyGlobalElements());
  mgr.setNodalColMap(numBoxNodes, nodalColMap.MyGlobalElements());
  mgr.finalizeBlockStructure();

  mgr.initializeProblemValues();

  for (int i=0; i<numDensity; i++) {
    int ownedPhysicsId = densityequ[i];
    for (int j=0; j<numOwnedNodes; j++) {
      int ownedNode = j;
      mgr.insertMatrixValue(ownedPhysicsId, ownedNode, ownedPhysicsID, ownedNode, 1.0); // insert diagonal
      for (int i1=0; i1<numUnknownsPerNode; i1++) {
	for (int j1=0; j1<numOwnedNodes; j1++) {
    

  mgr.insertRhsValue(i, j, rhsValue)

  // Construct a Map with NumElements and index base of 0
  Epetra_Map Map(NumElements, 0, Comm);


  // Create x and b vectors
  Epetra_Vector x(Map);
  Epetra_Vector b(Map);

  b.Random();
  x.Update(2.0, b, 0.0); // x = 2*b

  double bnorm, xnorm;
  x.Norm2(&xnorm);
  b.Norm2(&bnorm);

  cout << "2 norm of x = " << xnorm << endl
       << "2 norm of b = " << bnorm << endl;

  MPI_Finalize() ;

  return 0;
}

