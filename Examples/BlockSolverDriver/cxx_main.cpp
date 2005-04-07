/*@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
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
// ***********************************************************************
//@HEADER
*/

#include "Epetra_ConfigDefs.h"
#include "EpetraExt_Version.h"
#ifdef EPETRA_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_BlockMap.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"

#include "Epetra_CrsMatrix.h"
#include "EpetraExt_CrsMatrixIn.h"

int main(int argc, char *argv[]) {

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm comm (MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  int MyPID = comm.MyPID();

  bool verbose = false;
  bool verbose1 = false; 
  // Check if we should print results to standard out
  if (argc > 2) {
    if ((argv[2][0] == '-') && (argv[2][1] == 'v')) {
      verbose1 = true;
      if (MyPID==0) verbose = true;
    }
  }
  if (verbose)
    cout << EpetraExt::EpetraExt_Version() << endl << endl;

  if (verbose1) cout << comm << endl;


  // Uncomment the next three lines to debug in mpi mode
  //int tmp;
  //if (MyPID==0) cin >> tmp;
  //comm.Barrier();

  const int lineLength = 1025;
  char line[lineLength];
  int M, N, NZ;

  FILE * handle = 0;

  handle = fopen("A11","r");  // Open file

  // Strip off header lines (which start with "%")
  do {
    if(fgets(line, lineLength, handle)==0) return(-1);
  } while (line[0] == '%');

  // Get problem dimensions: M, N, NZ
  if(sscanf(line, "%d %d %d", &M, &N, &NZ)==0) return(-1);

  assert(M==N); // We need the matrix to be square.
  fclose(handle);

  Epetra_Map rowMap1(M, 0, comm); // Create default linear distribution

  handle = fopen("A22","r");  // Open file

  // Strip off header lines (which start with "%")
  do {
    if(fgets(line, lineLength, handle)==0) return(-1);
  } while (line[0] == '%');

  // Get problem dimensions: M, N, NZ
  if(sscanf(line, "%d %d %d", &M, &N, &NZ)==0) return(-1);

  assert(M==N); // We need the matrix to be square.
  fclose(handle);

  Epetra_Map rowMap2(M, 0, comm); // Create default linear distribution
  Epetra_CrsMatrix * A11; 
  Epetra_CrsMatrix * A12; 
  Epetra_CrsMatrix * A21; 
  Epetra_CrsMatrix * A22; 
  Epetra_Vector x1(rowMap1); x1.Random();
  Epetra_Vector x2(rowMap1); x2.Random();
  Epetra_Vector b1(rowMap1); b1.Random();
  Epetra_Vector b2(rowMap1); b2.Random();

  EpetraExt::MatrixMarketFileToCrsMatrix("A11", rowMap1, A11);
  EpetraExt::MatrixMarketFileToCrsMatrix("A12", rowMap1, A12);
  EpetraExt::MatrixMarketFileToCrsMatrix("A21", rowMap2, A21);
  EpetraExt::MatrixMarketFileToCrsMatrix("A22", rowMap2, A22);

  A11->OptimizeStorage();
  A12->OptimizeStorage();
  A21->OptimizeStorage();
  A22->OptimizeStorage();

  double residual;
  residual = A->NormInf(); double rAInf = residual;
  if (verbose) cout << "Inf Norm of A                                                    = " << residual << endl;
  residual = A->NormOne(); double rAOne = residual;
  if (verbose) cout << "One Norm of A                                                    = " << residual << endl;
  x.Norm2(&residual); double rAx = residual;
  if (verbose) cout << "Norm of A*x                                                     = " << residual << endl;
  b.Norm2(&residual); double rb = residual;
  if (verbose) cout << "Norm of b (should equal norm of Ax)                              = " << residual << endl;

  Epetra_Flops flopcounter;
  Epetra_Time timer(comm);

  A->SetFlopCounter(flopcounter);

  for (int j=0; j<2; j++) {
    bool TransA = (j==1);
    timer.ResetStartTime();

    if (!TransA)
      //10 matvecs
      for( int i = 0; i < 10; ++i )
	A->Multiply(TransA, x, b); // Compute b = A*x
    else
      //10 matvecs transpose
      for( int i = 0; i < 10; ++i )
	A->Multiply(TransA, b, x); // Compute x = A'*b
    
    double elapsed_time = timer.ElapsedTime();
    double total_flops = A->Flops();
    
    double MFLOPs = total_flops/elapsed_time/1000000.0;
    if (verbose) cout << "Total MFLOPs for 10 MatVec's with (Trans = " << TransA
		      << ") = " << MFLOPs << " (" << elapsed_time << " s)" <<endl;
    
  }
  
  delete A;


#ifdef EPETRA_MPI
  MPI_Finalize() ;
#endif

return 0 ;
}
