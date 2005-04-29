#include <mpi.h>
#include "dft_SolverManager.hpp"
#include <vector>
#include <cassert>
int main(int argc, char **argv)
{
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
  if (rank>0) {
    boxGIDs[ibox++] = ownedGIDs[0] - 2;
    boxGIDs[ibox++] = ownedGIDs[0] - 1;
  }
  for (int i=0; i<numOwned; i++) boxGIDs[ibox++] = ownedGIDs[i];
  if (rank<size-1) {
    boxGIDs[ibox++] = ownedGIDs[0] + 1;
    boxGIDs[ibox++] = ownedGIDs[0] + 2;
  }
  assert(ibox==numBox);

  mgr.setNodalRowMap(numOwned, &(ownedGIDs[0]));
  mgr.setNodalColMap(numBox, &(boxGIDs[0]));

  mgr.finalizeBlockStructure();

  mgr.initializeProblemValues();

  for (int j=0; j<numUnks; j++)
    for (int i=0; i<numOwned; i++) {
      mgr.insertRhsValue(j, i, (double) i*j);
      mgr.insertMatrixValue(j, i, j, i, -1.0);
    }

  mgr.finalizeProblemValues();
  
  mgr.setupSolver();

  mgr.solve();
  
  std::vector<double> allx(numUnks*numOwned);
  double ** x = new double*[numUnks];
  for (int i=0;i<numUnks; i++) x[i] = &(allx[0]) + i*numOwned;
  mgr.getLhs(x);

  //std::cout < allx << endl;

  MPI_Finalize();
  return(0);
}
