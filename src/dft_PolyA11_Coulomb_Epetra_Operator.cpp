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

#include "dft_PolyA11_Coulomb_Epetra_Operator.hpp"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_TestForException.hpp"

//=============================================================================
/*dft_PolyA11_Coulomb_Epetra_Operator::dft_PolyA11_Coulomb_Epetra_Operator(const Epetra_Map & ownedMap, const Epetra_Map & block1Map, const Epetra_Map & allGMap, const Epetra_Map & poissonMap, int * solverOptions, double * solverParams)
  : dft_PolyA11_Epetra_Operator(ownedMap, allGMap),
    solverOptions_(solverOptions),
    solverParams_(solverParams),
    allGMap_(allGMap),
    poissonMap_(poissonMap),
    block1Map_(block1Map),
    curPoissonRow_(-1) {
  Label_ = "dft_PolyA11_Coulomb_Epetra_Operator";
  poissonMatrix_ = new Epetra_CrsMatrix(Copy, poissonMap, 0);
  poissonMatrix_->SetLabel("PolyA11Coulomb::poissonMatrix");
  return;
  }*/
//=============================================================================
dft_PolyA11_Coulomb_Epetra_Operator::dft_PolyA11_Coulomb_Epetra_Operator(const Epetra_Map & ownedMap, const Epetra_Map & block1Map, const Epetra_Map & allGMap, const Epetra_Map & poissonMap, Teuchos::ParameterList * parameterList) 
  : dft_PolyA11_Epetra_Operator(ownedMap, allGMap),
//dft_PolyA11_Epetra_Operator(ownedMap, block1Map),
    parameterList_(parameterList),
    allGMap_(allGMap),
    poissonMap_(poissonMap),
    block1Map_(block1Map),
    curPoissonRow_(-1) {

  Label_ = "dft_PolyA11_Coulomb_Epetra_Operator";
  poissonMatrix_ = new Epetra_CrsMatrix(Copy, poissonMap, 0); //or ownedMap??
  poissonMatrix_->SetLabel("PolyA11Coulomb::poissonMatrix");
  return;
  }
//==============================================================================
dft_PolyA11_Coulomb_Epetra_Operator::~dft_PolyA11_Coulomb_Epetra_Operator() {
  /*  for (int i=0; i<numBlocks_; i++) if (matrix_[i]!=0) delete matrix_[i];
      delete [] matrix_;*/ //dft_PolyA11_Epetra_Operator destructor deletes these
  delete poissonMatrix_;
  return;
}
//=============================================================================
int dft_PolyA11_Coulomb_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    for (int i=0; i<numBlocks_; i++)
      matrix_[i]->PutScalar(0.0);
    poissonMatrix_->PutScalar(0.0);
  }  

  return(0);
}
//=============================================================================
int dft_PolyA11_Coulomb_Epetra_Operator::insertMatrixValue(int ownedPhysicsID, int ownedNode, int rowGID, int colGID, double value) {
  
  if (ownedPhysicsID >= numBlocks_) { //insert it into Poisson part
    if (firstTime_) {
      if (rowGID != curPoissonRow_) {
	insertPoissonRow();
	curPoissonRow_ = rowGID;
	curPoissonOwnedNode_ = ownedNode;
      }
      curPoissonRowValues_[colGID] += value;
    } 
    else {
      poissonMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
    }
  }
  else { //insert it into G part

  if (rowGID!=colGID) value = -value; // negate off-diagonal values to simplify kernel calls

  if (firstTime_) {
    if (rowGID!=curRow_) { 
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
      curOwnedPhysicsID_ = ownedPhysicsID;
      curOwnedNode_ = ownedNode;
    }
    curRowValues_[colGID] += value;
  }
  else
    matrix_[ownedPhysicsID]->SumIntoGlobalValues(ownedNode, 1, &value, &colGID);
  }
  return(0);
}
//=============================================================================
int dft_PolyA11_Coulomb_Epetra_Operator::insertPoissonRow() {
  if (curPoissonRowValues_.empty()) return(0);
  int numEntries = curPoissonRowValues_.size();
  if (numEntries>indices_.Length()) {
    indices_.Resize(numEntries);
    values_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos=curPoissonRowValues_.begin(); pos!=curPoissonRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  poissonMatrix_->InsertGlobalValues(curPoissonRow_, numEntries, values_.Values(), indices_.Values());
  curPoissonRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_PolyA11_Coulomb_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do
  
  if (firstTime_) { 
    insertRow(); // Dump any remaining entries
    insertPoissonRow();
  }
  for (int i=0; i<numBlocks_; i++) {
    matrix_[i]->FillComplete(allGMap_, ownedMap_);
    matrix_[i]->OptimizeStorage();
    //TEST_FOR_EXCEPT(!matrix_[i]->LowerTriangular());
  }
  poissonMatrix_->FillComplete(poissonMap_, poissonMap_); //or FillComplete();?
  poissonMatrix_->OptimizeStorage();

  isLinearProblemSet_ = true;
  firstTime_ = false;

  return(0);
}
//==============================================================================
int dft_PolyA11_Coulomb_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  int NumVectors = Y.NumVectors();
  int numMyElements = ownedMap_.NumMyElements();

  Y=X; // We can safely do this
  double ** Yptr;
  Y.ExtractView(&Yptr); //Get array of pointers to columns of Y

   double ** curY2 = new double *[NumVectors];
  for (int i = 0; i<NumVectors; i++) curY2[i] = Yptr[i]+numMyElements*numBlocks_;
  Epetra_MultiVector Y2(View, poissonMap_, curY2, NumVectors);

  Epetra_MultiVector Y1(View, allGMap_, Yptr, NumVectors);
  double ** curY1 = new double *[NumVectors];
    for (int i=0; i<NumVectors; i++) curY1[i] = Yptr[i];
    Epetra_MultiVector Y1tmp(View, ownedMap_, curY1, NumVectors); // Start Y1tmp to view first numNodes elements of Y1

  for (int i=0; i< numBlocks_; i++) {
        matrix_[i]->Multiply(false, Y1, Y1tmp);
	//matrix_[i]->Multiply(false, Y, Y1tmp);
    for (int j=0; j<NumVectors; j++) curY1[j]+=numMyElements; // Increment pointers to next block
    Y1tmp.ResetView(curY1); // Reset view to next block
  }

  //now to apply the inverse of poissonMatrix_ to Y2
  Epetra_LinearProblem implicitProblem(poissonMatrix_, &Y2, &Y2);

  int solverInt = Teuchos::getParameter<int>(*parameterList_, "Solver");
  //int solverInt = solverOptions_[AZ_solver];
  char * solverName;
  
  switch (solverInt) {
  case AM_lapack: solverName = "Amesos_Lapack"; break;
  case AM_klu: solverName = "Amesos_Klu"; break;
  case AM_mumps: solverName = "Amesos_Mumps"; break;
  case AM_umfpack: solverName = "Amesos_Umfpack"; break;
  case AM_superlu: solverName = "Amesos_Superlu"; break;
  case AM_superludist: solverName = "Amesos_Superludist"; break;
  case AM_pardiso: solverName = "Amesos_Pardiso"; break;
  case AM_taucs: solverName = "Amesos_Taucs"; break;
  default: break;
  }
  
  if (solverInt >= AM_lapack && solverInt <= AM_taucs) {
    EpetraExt::LinearProblem_Reindex::LinearProblem_Reindex(*ownedMap_);
    EpetraExt::LinearProblem_Reindex reindex(NULL);
    Epetra_LinearProblem reindexedProblem = reindex(implicitProblem);
    Amesos Amesos_Factory;
    Teuchos::RefCountPtr<Amesos_BaseSolver> directSolver_ = Teuchos::rcp(Amesos_Factory.Create(solverName, reindexedProblem));
    directSolver_->SetParameters(*parameterList_);
    directSolver_->SymbolicFactorization();
    directSolver_->NumericFactorization();
    directSolver_->Solve();
  } 
  else {
    Teuchos::RefCountPtr<AztecOO> solver_ = Teuchos::rcp(new AztecOO(implicitProblem));
    solver_->SetParameters(*parameterList_);
    //solver_->SetAllAztecOptions(solverOptions_);
    //solver_->SetAllAztecParams(solverParams_);
    solver_->Iterate(Teuchos::getParameter<int>(*parameterList_, "Max_iter"), Teuchos::getParameter<double>(*parameterList_, "Tol")); //Try to solve
    //solver_->Iterate(solverOptions_[AZ_max_iter], solverParams_[AZ_tol]);
    }

  delete [] curY1;
  delete [] curY2;
  return(0);
}
//==============================================================================
int dft_PolyA11_Coulomb_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numMyElements = ownedMap_.NumMyElements();

  double ** curY1 = new double *[NumVectors];
  double ** curY2 = new double *[NumVectors];
  double ** curX = new double *[NumVectors];
  double ** Yptr;
  double ** Xptr;

  Y.ExtractView(&Yptr); // Get array of pointers to columns of Y
  X.ExtractView(&Xptr); // Get array of pointers to columns of X

  Epetra_MultiVector X1(View, allGMap_, Xptr, NumVectors);
  Epetra_MultiVector X2(View, poissonMap_, Xptr, NumVectors);
  for (int i=0; i<NumVectors; i++) curY1[i] = Yptr[i];
  Epetra_MultiVector Y1tmp(View, ownedMap_, curY1, NumVectors); // Start Y1tmp to view first numNodes elements of Y1

  for (int i=0; i<NumVectors; i++) curX[i] = Xptr[i];
  Epetra_MultiVector Xtmp(View, ownedMap_, curX, NumVectors); // Start Xtmp to view first numNodes elements of X

  for (int i=0; i< numBlocks_; i++) {
    matrix_[i]->Multiply(false, X1, Y1tmp); // This gives a result that is X - off-diagonal-matrix*X
    Y1tmp.Update(-2.0, Xtmp, 1.0); // This gives a result of -X - off-diagonal-matrix*X
    Y1tmp.Scale(-1.0); // Finally negate to get the desired result
    for (int j=0; j<NumVectors; j++) {
      curY1[j]+=numMyElements; // Increment pointers to next block
      curX[j]+=numMyElements; // Increment pointers to next block
    }
    Y1tmp.ResetView(curY1); // Reset view to next block
    Xtmp.ResetView(curX); // Reset view to next block
  }

  //now to apply the poissonMatrix_ to the last chunk of X
  for (int i = 0; i<NumVectors; i++) curY2[i] = Yptr[i]+numMyElements*numBlocks_;
  Epetra_MultiVector Y2(View, poissonMap_, curY2, NumVectors);

  poissonMatrix_->Multiply(false, X2, Y2);

  delete [] curY1;
  delete [] curY2;
  delete [] curX;
  return(0);
}
