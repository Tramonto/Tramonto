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

//=============================================================================
dft_2x2_A11_invertible_solver::dft_2x2_A11_invertible_solver(int numBlocks, Epetra_CrsMatrix *** blockMatrix) 
  : numBlocks_(numBlocks),
    blockMatrix_(blockMatrix) {

  bool debug = true;

  if (debug) assert(blockMatrix-{i][j]->LowerTriangular());

  return;
}

//=============================================================================
dft_SolverManager::~dft_SolverManager() {
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
