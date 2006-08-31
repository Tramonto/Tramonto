//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
//@HEADER

#include <mpi.h>
#include "dft_SolverManager.hpp"
#include <vector>
#include <cassert>
int main(int argc, char **argv)
{

  bool debug = true;

  MPI_Init(&argc,&argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int numUnks = 3;
  dft_SolverManager mgr(numUnks, 0, 0, 0, MPI_COMM_WORLD);

  int numOwned = 4; 
  assert(numOwned>3);
  int numBox = numOwned;
  if (rank>0) numBox +=2;
  if (rank<size-1) numBox +=2;
  std::vector<int> ownedGIDs(numOwned);
  for (int i=0; i<numOwned; i++) ownedGIDs[i] = numOwned*rank+i;

  std::vector<int> boxGIDs(numBox);
  int ibox = 0;
  int ownedOffset = 0;
  if (rank>0) {
    boxGIDs[ibox++] = ownedGIDs[0] - 2;
    boxGIDs[ibox++] = ownedGIDs[0] - 1;
    ownedOffset = 2;
  }
  for (int i=0; i<numOwned; i++) boxGIDs[ibox++] = ownedGIDs[i];
  if (rank<size-1) {
    boxGIDs[ibox++] = ownedGIDs[numOwned-1] + 1;
    boxGIDs[ibox++] = ownedGIDs[numOwned-1] + 2;
  }
  assert(ibox==numBox);

  if (debug) std::cout << "Calling setNodalRowMap" <<std::endl;
  mgr.setNodalRowMap(numOwned, &(ownedGIDs[0]));

  if (debug) std::cout << "Calling setNodalColMap" <<std::endl;
  mgr.setNodalColMap(numBox, &(boxGIDs[0]));

  if (debug) std::cout << "Calling finalizeBlockStructure" <<std::endl;
  mgr.finalizeBlockStructure();

  if (debug) std::cout << "Calling initializeProblemValues" <<std::endl;
  mgr.initializeProblemValues();

  if (debug) std::cout << "Calling insertRhsValue and insertMatriValue" <<std::endl;
  for (int j=0; j<numUnks; j++)
    for (int i=0; i<numOwned; i++) {
      mgr.insertRhsValue(j, i, (double) 1+i+numOwned*j +numUnks*numOwned*rank);
      mgr.insertMatrixValue(j, i, j, i+ownedOffset, -1.0);
    }

  if (debug) std::cout << "Calling finalizeProblemValues" <<std::endl;
  mgr.finalizeProblemValues();
  
  std::vector<double> allb(numUnks*numOwned);
  double ** b = new double*[numUnks];
  for (int i=0;i<numUnks; i++) b[i] = &(allb[0]) + i*numOwned;
  if (debug) std::cout << "Calling getRhs" <<std::endl;
  mgr.getRhs(b);

  for (int i=0; i<numUnks*numOwned; i++) std::cout << allb[i] << endl;

  if (debug) std::cout << "Calling setupSolver" <<std::endl;
  mgr.setupSolver();

  if (debug) std::cout << "Calling solve" <<std::endl;
  mgr.solve();
  
  std::vector<double> allx(numUnks*numBox);
  double ** x = new double*[numUnks];
  for (int i=0;i<numUnks; i++) x[i] = &(allx[0]) + i*numBox;
  if (debug) std::cout << "Calling getLhs" <<std::endl;
  mgr.getLhs(x);

  for (int i=0; i<numUnks; i++) 
    for (int j=0; j<numBox; j++) 
      std::cout << "Rank:"<< rank << " Unk:"<< i << " BoxID:"<<j<< " x:" << x[i][j] << endl;

  MPI_Finalize();
  return(0);
}
