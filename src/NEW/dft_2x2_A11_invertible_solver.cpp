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
  : op_(blockMatrix[1][1], blockMatrix[1][2], blockMatrix[2][1], blockMatrix[2][2]),
    numBlocks_(numBlocks),
    blockMatrix_(blockMatrix) {

  bool debug = true;

  if (debug) assert(blockMatrix-{i][j]->LowerTriangular());

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
