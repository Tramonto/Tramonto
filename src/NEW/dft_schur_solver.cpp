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

#include "mpi.h"
#include "Epetra_Comm.h"
#include "Epetra_Object.h"
#include "Epetra_MpiComm.h"
#include "az_aztec.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h"
#include "AztecOO.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_Time.h"
#include "dft_schur_epetra_operator.h"
#include "dft_globals_const.h"
#include "rf_allo.h"

#include "dft_schur_solver.h"
#include <cassert>

void dft_create_schur_solver(int * proc_config, int * external, int * bindx, 
			     double * val, int * update, int * update_index,
			     int * extern_index, int * data_org, int N_update, 
			     int Nunk_per_node,
			     DFT_SCHUR_SOLVER ** schur_solver) {

  Epetra_MpiComm comm(MPI_COMM_WORLD);
  Epetra_Time timer(comm);
  bool debug = false;

  int nnodes = N_update/Nunk_per_node;
  assert(nnodes*Nunk_per_node==N_update); // Make sure there is a constant number of unknowns per node

  int nexternal = data_org[AZ_N_external];
  int * ColGIDs = new int[N_update+nexternal];
  if (debug) for (int i=0; i<N_update+nexternal; i++) ColGIDs[i] = -1;
  for (int i=0; i<N_update; i++) ColGIDs[update_index[i]] = update[i];
  for (int i=0; i<nexternal; i++) ColGIDs[extern_index[i]] = external[i];
  if (debug) for (int i=0; i<N_update+nexternal; i++) assert(ColGIDs[i]!=-1); // Make sure all ColGIDs were assigned
  if (debug) for (int i=0; i<N_update; i++) assert(update_index[i]==i); // We assume that local IDs are ordered 0 to N_update-1
  if (debug) for (int i=0; i<nexternal; i++) assert(extern_index[i]==N_update+i);

  int * node_reorder = new int[Nunk_per_node];
  
  /* These parameter values are for the polymer problems
  int n2 = 6;
  int n1 = Nunk_per_node-n2;
  int nforward = n1/2;
  for (int i=0; i<nforward; i++) {
    node_reorder[i] = n2 + i*2;
    node_reorder[i] = Nunk_per_node - i*2;
  }

  for (int i=0; i<n2; i++) node_reorder[n1+i] = i;
  */

  // These parameters are for a hard sphere with implicit nonlocal densities.
  int n2 = 1;
  int n1 = Nunk_per_node-n2;
  node_reorder[0] = 1;
  node_reorder[1] = 2;
  node_reorder[2] = 5;
  node_reorder[3] = 6;
  node_reorder[4] = 7;
  node_reorder[5] = 3;
  node_reorder[6] = 4;
  node_reorder[7] = 8;
  node_reorder[8] = 9;
  node_reorder[9] = 10;
  node_reorder[10] = 0;

  //Common code, using node_order as defined above
  
  int nrows1 = n1*nnodes;
  int nrows2 = N_update - nrows1;


  int * reorder = new int[N_update];
  int * reorder_cols1 = new int[nrows1+nexternal];
  int * reorder_cols2 = new int[nrows2+nexternal];

  for (int j=0; j<nnodes; j++) {
    int jinc = update[j*Nunk_per_node];
    for (int i=0; i<Nunk_per_node; i++) {
      reorder[nnodes*i+j] = jinc + node_reorder[i];
      assert(update[j*Nunk_per_node+i]==update[j*Nunk_per_node]+i); // Check to see that GIDs at a node are contiguously ordered
    }
  }

  for (int i=0; i<nrows1; i++) reorder_cols1[i] = reorder[i];
  for (int i=0; i<nrows2; i++) reorder_cols2[i] = reorder[nrows1+i];
  int nexternal1 = 0;
  int nexternal2 = 0;
  for (int i=0; i<nexternal; i++) {
    int nodeID = ColGIDs[extern_index[i]];
    if (nodeID%Nunk_per_node<n2) { // Type 2 column
      reorder_cols2[nrows2+nexternal2] = nodeID;
      nexternal2++;
    }
    else {
      reorder_cols1[nrows1+nexternal1] = nodeID;
      nexternal1++;
    }
  }
  int ncols1 = nrows1 + nexternal1;
  int ncols2 = nrows2 + nexternal2;
      

  Epetra_Map  RowMap1(-1, nrows1, reorder, 0, comm);
  Epetra_Map  RowMap2(-1, nrows2, reorder+nrows1, 0, comm);
  Epetra_Map  ColMap1(-1, ncols1, reorder_cols1, 0, comm);
  Epetra_Map  ColMap2(-1, ncols2, reorder_cols2, 0, comm);
  Epetra_CrsMatrix * A11 = new Epetra_CrsMatrix(Copy, RowMap1, ColMap1, 0);
  Epetra_CrsMatrix * A12 = new Epetra_CrsMatrix(Copy, RowMap1, ColMap2, 0);
  Epetra_CrsMatrix * A21 = new Epetra_CrsMatrix(Copy, RowMap2, ColMap1, 0);
  Epetra_CrsMatrix * A22 = new Epetra_CrsMatrix(Copy, RowMap2, ColMap2, 0);

  Epetra_Map  RowMap(-1, N_update, update, 0, comm);
  Epetra_Map  ColMap(-1, N_update+nexternal, ColGIDs, 0, comm);
  Epetra_CrsMatrix * A = new Epetra_CrsMatrix(Copy, RowMap, ColMap, 0);

  int * curColGIDs = new int[N_update+nexternal+1];
  double * curGlobalVals = new double[N_update+nexternal+1];

  for (int i=0; i< N_update; i++) {
    int curRowGID = update[i];
    int curNumEntries = bindx[i+1] - bindx[i];
    int * curEntries = bindx+bindx[i];
    double * curVals = val+bindx[i];

    for (int j=0; j<curNumEntries; j++) curGlobalVals[j] = curVals[j];
    curGlobalVals[curNumEntries] = val[i]; // Insert diagonal value

    for (int j=0; j<curNumEntries; j++) curColGIDs[j] = ColGIDs[curEntries[j]];
    curColGIDs[curNumEntries] = curRowGID; // Insert diagonal index

    if (RowMap1.MyGID(curRowGID)) {
      A11->InsertGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
      A12->InsertGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
    }
    else {
      A21->InsertGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
      A22->InsertGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
    }
    A->InsertGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
  }

  A11->FillComplete();
  A22->FillComplete();
  A->FillComplete();
  A12->FillComplete(RowMap2, RowMap1);
  A21->FillComplete(RowMap1, RowMap2);
  A11->OptimizeStorage();
  A12->OptimizeStorage();
  A21->OptimizeStorage();
  A22->OptimizeStorage();

  if (debug) assert(A11->LowerTriangular());

  *schur_solver = (DFT_SCHUR_SOLVER *) new DFT_SCHUR_SOLVER;
  (*schur_solver)->A11 = (void *) A11;
  (*schur_solver)->A12 = (void *) A12;
  (*schur_solver)->A21 = (void *) A21;
  (*schur_solver)->A22 = (void *) A22;
  (*schur_solver)->A   = (void *) A;
  
  delete [] ColGIDs;
  delete [] node_reorder;
  delete [] reorder;
  delete [] reorder_cols1;
  delete [] reorder_cols2;
  delete [] curColGIDs;
  delete [] curGlobalVals;

  cout << "Time in dft_create_schur_solver = " << timer.ElapsedTime();
  return;
}

void dft_update_schur_solver(int * bindx, double * val,
			     DFT_SCHUR_SOLVER * schur_solver) {

  Epetra_CrsMatrix * A11 = (Epetra_CrsMatrix *) schur_solver->A11;
  Epetra_CrsMatrix * A12 = (Epetra_CrsMatrix *) schur_solver->A12;
  Epetra_CrsMatrix * A21 = (Epetra_CrsMatrix *) schur_solver->A21;
  Epetra_CrsMatrix * A22 = (Epetra_CrsMatrix *) schur_solver->A22;
  Epetra_CrsMatrix * A   = (Epetra_CrsMatrix *) schur_solver->A;

  Epetra_Time timer(A11->Comm());

  int * update = A->RowMap().MyGlobalElements();
  int * ColGIDs = A->ColMap().MyGlobalElements();
  int N_update = A->RowMap().NumMyElements();
  int N_cols = A->ColMap().NumMyElements();

  Epetra_Map RowMap1(A11->RowMap());

  int * curColGIDs = new int[N_cols];
  double * curGlobalVals = new double[N_cols];

  for (int i=0; i< N_update; i++) {
    int curRowGID = update[i];
    int curNumEntries = bindx[i+1] - bindx[i];
    int * curEntries = bindx+bindx[i];
    double * curVals = val+bindx[i];

    for (int j=0; j<curNumEntries; j++) curGlobalVals[j] = curVals[j];
    curGlobalVals[curNumEntries] = val[i]; // Insert diagonal value

    for (int j=0; j<curNumEntries; j++) curColGIDs[j] = ColGIDs[curEntries[j]];
    curColGIDs[curNumEntries] = curRowGID; // Insert diagonal index

    if (RowMap1.MyGID(curRowGID)) {
      A11->ReplaceGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
      A12->ReplaceGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
    }
    else {
      A21->ReplaceGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
      A22->ReplaceGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
    }
    A->ReplaceGlobalValues(curRowGID, curNumEntries+1, curGlobalVals, curColGIDs);
  }
  delete [] curColGIDs;
  delete [] curGlobalVals;

  cout << "Time in dft_update_schur_solver = " << timer.ElapsedTime();
  return;
}

void dft_apply_schur_solver(DFT_SCHUR_SOLVER * schur_solver, double * x_aztec, double * b_aztec) {

  bool debug = false;

  Epetra_CrsMatrix * A11 = (Epetra_CrsMatrix *) schur_solver->A11;
  Epetra_CrsMatrix * A12 = (Epetra_CrsMatrix *) schur_solver->A12;
  Epetra_CrsMatrix * A21 = (Epetra_CrsMatrix *) schur_solver->A21;
  Epetra_CrsMatrix * A22 = (Epetra_CrsMatrix *) schur_solver->A22;
  Epetra_CrsMatrix * A   = (Epetra_CrsMatrix *) schur_solver->A;

  Epetra_Time timer(A11->Comm());

  if (debug) {
    EpetraExt::RowMatrixToMatrixMarketFile("A11.mm", *A11, "A11 Matrix",
					   "A11 matrix from Tramonto Schur complement");
    EpetraExt::RowMatrixToMatrixMarketFile("A12.mm", *A12, "A12 Matrix",
					   "A12 matrix from Tramonto Schur complement");
    EpetraExt::RowMatrixToMatrixMarketFile("A21.mm", *A21, "A21 Matrix",
					   "A21 matrix from Tramonto Schur complement");
    EpetraExt::RowMatrixToMatrixMarketFile("A22.mm", *A22, "A22 Matrix",
					   "A22 matrix from Tramonto Schur complement");
    EpetraExt::RowMatrixToMatrixMarketFile("A.mm", *A, "A Matrix",
					   "A matrix from Tramonto Schur complement");
  }
  Epetra_Vector x(View, A->DomainMap(), x_aztec);
  Epetra_Vector b(View, A->RangeMap(), b_aztec);

  if (debug) {
    Epetra_Vector resid(b.Map());
    A->Multiply(false, x, resid);
    resid.Update(1.0, b, -1.0);
    double resval;
    resid.Norm2(&resval);
    cout << "Initial residual before Schur Solver = "<<resval<< endl;
  }

  Epetra_Map RowMap(A->RowMap());
  Epetra_Map RowMap1(A11->RowMap());
  Epetra_Map RowMap2(A22->RowMap());
  
  Epetra_Vector x1(RowMap1);
  Epetra_Vector b1(RowMap1);
  Epetra_Vector x2(RowMap2);
  Epetra_Vector b2(RowMap2);
  Epetra_Vector b2s(RowMap2);

  Epetra_Import importer1(RowMap1, RowMap);
  Epetra_Import importer2(RowMap2, RowMap);

  x1.Import(x, importer1, Insert);
  b1.Import(b, importer1, Insert);
  x2.Import(x, importer2, Insert);
  b2.Import(b, importer2, Insert);
  
  dft_schur_epetra_operator op(A11, A12, A21, A22);
  op.ComputeRHS(b1, b2, b2s);
  Epetra_LinearProblem problem(&op, &x2, &b2s);
  AztecOO solver(problem);
  solver.SetPrecMatrix(A22);
  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  solver.Iterate(200, 1.0E-7);

  op.ComputeX1(b1, x2, x1);
  
  x.Export(x1, importer1, Insert);
  x.Export(x2, importer2, Insert);

  if (debug) {
    Epetra_Vector resid(b.Map());
    A->Multiply(false, x, resid);
    resid.Update(1.0, b, -1.0);
    double resval;
    resid.Norm2(&resval);
    cout << "Residual from Schur Solver = "<<resval<< endl;
  }
  cout << "Time in dft_update_schur_solver = " << timer.ElapsedTime();
}
void dft_destroy_schur_solver(DFT_SCHUR_SOLVER ** schur_solver) {

  Epetra_CrsMatrix * A11 = (Epetra_CrsMatrix *) (*schur_solver)->A11;
  Epetra_CrsMatrix * A12 = (Epetra_CrsMatrix *) (*schur_solver)->A12;
  Epetra_CrsMatrix * A21 = (Epetra_CrsMatrix *) (*schur_solver)->A21;
  Epetra_CrsMatrix * A22 = (Epetra_CrsMatrix *) (*schur_solver)->A22;
  Epetra_CrsMatrix * A   = (Epetra_CrsMatrix *) (*schur_solver)->A;

  delete A11;
  delete A12;
  delete A21;
  delete A22;
  delete A;

  delete (DFT_SCHUR_SOLVER *) (*schur_solver);
  *schur_solver = 0;
  return;
}
