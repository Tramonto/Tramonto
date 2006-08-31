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

#include "dft_2x2_A11_invertible_solver.h"
#include <cassert>

//=============================================================================
dft_2x2_A11_invertible_solver::dft_2x2_A11_invertible_solver(int numBlocks, Epetra_Operator *** blockMatrix) 
  : op_(blockMatrix[0][0], blockMatrix[0][1], blockMatrix[1][0], blockMatrix[1][1]),
    numBlocks_(numBlocks),
    blockMatrix_(blockMatrix) {

  bool debug = true;

  if (debug) assert(blockMatrix->[0][0]->LowerTriangular());

  return;
}

//=============================================================================
dft_2x2_A11_invertible_solver::~dft_SolverManager() {
  return;
}
//=============================================================================
dft_2x2_A11_invertible_solver::solve(Epetra_Vector ** blockRhs, Epetra_Vector ** blockLhs) {

  bool debug = true;
  bool writefiles = false;

  if (writefiles) {
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
  
  if (debug) {
    double residual = 0.0;
    residual = op_.ComputeGlobalResidual(*blockRhs[1], *blockRhs[2], *blockLhs[1], *blockLhs[2]);
    if (blockLhs[1]->Comm().MyPID()==0) cout << "Initial residual = " << residual << endl;
  }
  
  Epetra_Vector b2s(blockRhs[2]->Map());
  op_.ComputeRHS(*blockRhs[1], *blockRhs[2], b2s);
  Epetra_LinearProblem problem(&op, blockLhs[2], &b2s);
  AztecOO solver(problem);
  solver.SetPrecMatrix(blockMatrix_[2][2]);
  solver.SetAztecOption(AZ_precond, AZ_Jacobi);
  solver.Iterate(200, 1.0E-7);

  op.ComputeX1(*blockRhs[1], *blockLhs[2], *blockLhs[1]);
  
  if (debug) {
    double residual = 0.0;
    residual = op_.ComputeGlobalResidual(*blockLhs[1], *blockLhs[2], *blockRhs[1], *blockRhs[2]);
    if (blockLhs[1]->Comm().MyPID()==0) cout << "Final residual = " << residual << endl;
  }
}
