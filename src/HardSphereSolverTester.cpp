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

#include "dft_HardSphereLinProbMgr.hpp"
#include <mpi.h>

//int gid(int physicsID, nodeID, numPhysics, numNodes) {return(
int main(int argc, char *argv[])
{


  // Initialize MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm( MPI_COMM_WORLD );

  //if (comm.NumProc()>1) { cout << "Can only use a single processor for now!!" << endl; exit(1);}
  //int tmp;
  //if (comm.MyPID()==0) std::cin >> tmp;
  //comm.Barrier();

  int numOwnedNodes = 4;

  
  
  int numUnknownsPerNode = 5;
  int densityequ[] = {0, 1}; int numDensity = 2;
  int indnonlocalequ[] = {2, 3, 4}; int numIndNonLocal = 3;
  int depnonlocalequ[] = {5, 6, 7}; int numDepNonLocal = 0;
  //int depnonlocalequ[] = {2, 4, 6}; int numDepNonLocal = 3;
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
  

  // For now we have a trivial box
  int numBoxNodes = numOwnedNodes;
  Epetra_Map nodalColMap(nodalRowMap);

  dft_BasicLinProbMgr * basicMgr = new dft_BasicLinProbMgr(numUnknownsPerNode, 0, 0, MPI_COMM_WORLD);
  dft_HardSphereLinProbMgr * hardSphereMgr = new dft_HardSphereLinProbMgr(numUnknownsPerNode, 0, 0, MPI_COMM_WORLD, true);  // debug=true, build globalmatrix also.

  Epetra_SerialDenseVector densityDiagonalValues(numDensity); densityDiagonalValues.Random();
  Epetra_SerialDenseVector nonLocalOffDiagonalValues(numDepNonLocal); nonLocalOffDiagonalValues.Random();
  Epetra_SerialDenseVector densityOnNonLocal(numIndNonLocal+numDepNonLocal); densityOnNonLocal.Random();
  for (int imgr = 0; imgr < 2; imgr++) {
    dft_BasicLinProbMgr * mgr;
    if (imgr==0)
      mgr = basicMgr;
    else
      mgr = dynamic_cast<dft_BasicLinProbMgr *>(hardSphereMgr);
    assert(mgr!=0); // Make sure we got a real object

    mgr->setNodalRowMap(numOwnedNodes, nodalRowMap.MyGlobalElements());
    mgr->setNodalColMap(numBoxNodes, nodalColMap.MyGlobalElements());

    if (imgr==1) {  // These are specific to the HardSphereLinProbMgr
      hardSphereMgr->setIndNonLocalEquationIDs(numIndNonLocal, indnonlocalequ);
      hardSphereMgr->setDepNonLocalEquationIDs(numDepNonLocal, depnonlocalequ);
      hardSphereMgr->setDensityEquationIDs(numDensity, densityequ);
    }
    mgr->finalizeBlockStructure();
  
    for (int iters=0; iters<2; iters++) {

      mgr->initializeProblemValues();
      int irhs = 0;
      for (int i=0; i<numOwnedNodes; i++) {
	int ownedNode = i;

	// Density Equations
	for (int j=0; j<numDensity; j++) {
	  int ownedPhysicsID = densityequ[j];
	  mgr->insertRhsValue(ownedPhysicsID, ownedNode, (double) irhs++); // rhs
	  int boxPhysicsID = ownedPhysicsID; // density on density
	  int boxNode = ownedNode;
	  mgr->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, densityDiagonalValues[j]);
	  assert(mgr->getMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode)==densityDiagonalValues[j]);
      

	  int ratioNonLocaltoDensity = (numIndNonLocal+numDepNonLocal)/numDensity;
	  for (int k=0; k<ratioNonLocaltoDensity; k++) {
	    int k1 = j*ratioNonLocaltoDensity+k;
	    if (k1<numIndNonLocal)
	      boxPhysicsID = indnonlocalequ[k1];  // density on independent nonlocal density
	    else
	      boxPhysicsID = depnonlocalequ[k1-numIndNonLocal];  // density on dependent nonlocal density
	  
	    mgr->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
	    assert(mgr->getMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode)==1.0);
	  } // End of Density Equations
	}
    
	// Independent nonlocal Equations

	for (int j=0; j<numIndNonLocal; j++) {
	  int ownedPhysicsID = indnonlocalequ[j];
	  mgr->insertRhsValue(ownedPhysicsID, ownedNode, (double) irhs++); // rhs
	  int boxPhysicsID = ownedPhysicsID; // ind nonlocal on self
	  int boxNode = ownedNode;
	  mgr->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
	  assert(mgr->getMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode)==1.0);

	  boxPhysicsID = densityequ[j*numDensity/numIndNonLocal];  // indnonlocal on density
	  double value = 1 - 0.5*((double) abs(j))/((double) numIndNonLocal);
	  mgr->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
	  assert(mgr->getMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode)==value);
	} // End of independent nonlocal Equations
      
	// Denpendent nonlocal Equations
    
	for (int j=0; j<numDepNonLocal; j++) {
	  int ownedPhysicsID = depnonlocalequ[j];
	  mgr->insertRhsValue(ownedPhysicsID, ownedNode, (double) irhs++); // rhs
	  int boxPhysicsID = ownedPhysicsID; // dep nonlocal on self
	  int boxNode = ownedNode;
	  mgr->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, 1.0);
	  assert(mgr->getMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode)==1.0);

	  if (j>0) { // block subdiagonal 
	    boxPhysicsID = indnonlocalequ[j];  // dep on ind nonlocal
	    mgr->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, -0.5);
	  }
	  boxPhysicsID = densityequ[numDensity/numDepNonLocal]; // g on cms
	  double value = 1 - 0.5*((double) abs(j))/((double) numDepNonLocal);
	  mgr->insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode, value);
	  assert(mgr->getMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNode)==value);
	
	} // End of dependent nonlocal Equations
      } // End of owned nodes loop
      mgr->finalizeProblemValues();  // All done filling values
    } // End of iters loop
    
    mgr->writeMatrix("globalMatrix", "Matrix tester global matrix", "MatrixTesterGlobalMatrix");
    double * xptr = new double[numUnknownsPerNode*numBoxNodes];
    for (int i=0; i<numUnknownsPerNode*numBoxNodes; i++) xptr[i] = (double) i+1 ;
    double * bptr = new double[numUnknownsPerNode*numOwnedNodes];
    for (int i=0; i<numUnknownsPerNode*numOwnedNodes; i++) bptr[i] = 0.0;
    double ** x = new double *[numUnknownsPerNode];
    double ** b = new double *[numUnknownsPerNode];
    for (int i=0; i<numUnknownsPerNode; i++) x[i] = xptr+i*numBoxNodes;
    for (int i=0; i<numUnknownsPerNode; i++) b[i] = bptr+i*numOwnedNodes;
    
    if (imgr==1) // Peculiar to HardSphereLinProbMgr
      hardSphereMgr->Check(true);
    //mgr->applyMatrix((const double **) x, b);
    mgr->getRhs(b);
    double * tmp = bptr;
    for (int i=0; i<numUnknownsPerNode; i++)
      for (int j=0; j<numOwnedNodes; j++) std::cout << "b[physics="<<i<<"][node="<<nodalRowMap.GID(j)<<"] = " << *tmp++ << std::endl;
    
    mgr->setupSolver();
    for (int i=0; i<numUnknownsPerNode*numBoxNodes; i++) xptr[i] = 0.0 ;

    mgr->solve();

    mgr->getLhs(x);
    tmp = xptr;
    for (int i=0; i<numUnknownsPerNode; i++)
      for (int j=0; j<numBoxNodes; j++) std::cout << "x[physics="<<i<<"][node="<<nodalColMap.GID(j)<<"] = " << *tmp++ << std::endl;

    delete[] x;
    delete [] xptr;
    delete[] b;
    delete [] bptr;

  }// End imgr loop

  MPI_Finalize();

  return 0;
}

