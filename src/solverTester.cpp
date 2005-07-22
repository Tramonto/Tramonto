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
  int cmsIntegrationRange = 5; // Number of grid points involved in "F" block stencil.

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

  Epetra_Map globalRowMap(-1, numMyUnknowns, globalGids.Values(), 0, comm);
  Epetra_CrsMatrix(Copy, globalRowMap, 0);
  

  // For now we are on a single processor
  int numBoxNodes = numOwnedNodes;
  Epetra_Map nodalColMap(nodalRowMap);

  dft_PolyLinProbMgr mgr(numUnknownsPerNode, 0, 0, MPI_COMM_WORLD);
  mgr.setNodalRowMap(numOwnedNodes, nodalRowMap.MyGlobalElements());
  mgr.setNodalColMap(numBoxNodes, nodalColMap.MyGlobalElements());

  mgr.setGEquationIDs(numG, gequ);
  mgr.setGInvEquationIDs(numGinv, ginvequ);
  mgr.setCmsEquationIDs(numCms, cmsequ);
  mgr.setDensityEquationIDs(numDensity, densityequ);
  mgr.finalizeBlockStructure();

  mgr.initializeProblemValues();
  assert(numG==numGinv && numDensity==numCms);  // Sanity test for assumptions below
  for (int i=0; i<numOwnedNodes; i++) {
    int ownedNode = nodalRowMap.GID(i);

    // Density Equations

    for (int j=0; j<numDensity; j++) {
      int ownedPhysicsID = densityequ[j];
      mgr.insertRhsValue(ownedPhysicsID, ownedNode, 1.0); // rhs
      int boxPhysicsID = ownedPhysicsID; // density on density
      int boxNode = ownedNode;
      mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
      
      boxPhysicsID = cmsequ[j];  // density on cms
      mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0e-11);
      
      int ratioGtoDensity = numGinv/numDensity;
      for (int k=0; k<ratioGtoDensity; k++) {
	int k1 = j*ratioGtoDensity+k;
	boxPhysicsID = gequ[k1];  // density on G
	mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
	boxPhysicsID = ginvequ[k1];  // density on GInv
	mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
      } // End of Density Equations
    }
    
    // Cms Field Equations

    for (int j=0; j<numCms; j++) {
      int ownedPhysicsID = cmsequ[j];
      mgr.insertRhsValue(ownedPhysicsID, ownedNode, 1.0); // rhs
      int boxPhysicsID = ownedPhysicsID; // cms on cms
      int boxNode = ownedNode;
      mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
      
      for (int k=-cmsIntegrationRange; k<=cmsIntegrationRange; k++) {
	boxNode = ownedNode + k;
	boxPhysicsID = densityequ[j];  // cms on density
	if (boxNode>=0 && boxNode<numOwnedNodes) {
	  double value = 1 - 0.5*((double) abs(k))/((double) cmsIntegrationRange); // 1 on diagonal, taper off away from diagonal
	  mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
	}
      }
    } // End of CMS Equations
      
    // G Equations
    
    for (int j=0; j<numG; j++) {
      int ownedPhysicsID = gequ[j];
      mgr.insertRhsValue(ownedPhysicsID, ownedNode, 1.0); // rhs
	int boxPhysicsID = ownedPhysicsID; // g on g
      int boxNode = ownedNode;
      mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
      
      if (j>0) { // block subdiagonal 
	boxPhysicsID = gequ[j-1];  // g_j on g_{j-1}
	mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, -0.5);
      }
      boxPhysicsID = cmsequ[j/numCms]; // g on cms
      mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, -0.5);

    } // End of G Equations

    // G Inverse Equations

    for (int j=0; j<numGinv; j++) {
      int ownedPhysicsID = ginvequ[j];
      mgr.insertRhsValue(ownedPhysicsID, ownedNode, 1.0); // rhs
      int boxPhysicsID = ownedPhysicsID; // ginv on ginv
      int boxNode = ownedNode;
      mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
      
      if (j>0) { // block subdiagonal 
	boxPhysicsID = ginvequ[j-1];  // ginv_j on ginv_{j-1}
	mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, -0.5);
      }
      boxPhysicsID = gequ[j]; // ginv on g
      mgr.insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, -0.5);

    } // End of G Inverse Equations

  } // End of numOwnedNodes loop

  mgr.finalizeProblemValues();  // All done filling values

  mgr.setupSolver();

  mgr.solve();

  double * xptr = new double[numUnknownsPerNode*numOwnedNodes];
  double ** x = new double *[numUnknownsPerNode];
  for (int i=0; i<numUnknownsPerNode; i++) x[i] = xptr+i*numOwnedNodes;

  mgr.getLhs(x);

  delete[] x;
  delete [] xptr;

  MPI_Finalize();

  return 0;
}

