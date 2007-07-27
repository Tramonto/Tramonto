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

#include "dft_PolyA22_Coulomb_Epetra_Operator.hpp"
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

//============================================================================
/*dft_PolyA22_Coulomb_Epetra_Operator::dft_PolyA22_Coulomb_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & poissonMap, const Epetra_Map & cmsDensMap, const Epetra_Map & block2Map, int * options, double * params) 
  : dft_PolyA22_Epetra_Operator(cmsMap, densityMap, cmsDensMap, options, params),
    poissonMap_(poissonMap),
    cmsDensMap_(cmsDensMap),
    block2Map_(block2Map)
 {
  poissonMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, poissonMap, 0));
  cmsOnPoissonMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cmsMap, 0));
  poissonOnDensityMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, poissonMap, 0));
  Label_ = "dft_PolyA22_Coulomb_Epetra_Operator";
  cmsOnDensityMatrix_->SetLabel("PolyA22Coulomb::cmsOnDensityMatrix");
  poissonMatrix_->SetLabel("PolyA22Coulomb::poissonMatrix");
  cmsOnPoissonMatrix_->SetLabel("PolyA22Coulomb::cmsOnPoissonMatrix");
  poissonOnDensityMatrix_->SetLabel("PolyA22Coulomb::poissonOnDensityMatrix");
  ML_Epetra::SetDefaults("SA",MLList_);
  MLList_.set("ML output", 0);
  }*/
//=============================================================================
dft_PolyA22_Coulomb_Epetra_Operator::dft_PolyA22_Coulomb_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & poissonMap, const Epetra_Map & cmsDensMap, const Epetra_Map & block2Map, Teuchos::ParameterList * parameterList) 
  : dft_PolyA22_Epetra_Operator(cmsMap, densityMap, cmsDensMap, parameterList),
    poissonMap_(poissonMap),
    cmsDensMap_(cmsDensMap),
    block2Map_(block2Map),
    curPoissonRow_(-1),
    curCPRow_(-1),
    curPDRow_(-1)
 {

    poissonMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, poissonMap, 0));
    cmsOnPoissonMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cmsMap, 0));
    poissonOnDensityMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, poissonMap, 0));
   Label_ = "dft_PolyA22_Coulomb_Epetra_Operator";
   cmsOnDensityMatrix_->SetLabel("PolyA22Coulomb::cmsOnDensityMatrix");
   poissonMatrix_->SetLabel("PolyA22Coulomb::poissonMatrix");
   cmsOnPoissonMatrix_->SetLabel("PolyA22Coulomb::cmsOnPoissonMatrix");
   poissonOnDensityMatrix_->SetLabel("PolyA22Coulomb::poissonOnDensityMatrix");
   ML_Epetra::SetDefaults("SA",MLList_);
   MLList_.set("ML output", 0);
}
//==============================================================================
dft_PolyA22_Coulomb_Epetra_Operator::~dft_PolyA22_Coulomb_Epetra_Operator() {
  return;
}
//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::initializeProblemValues() {

  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    if (!isFLinear_) cmsOnDensityMatrix_->PutScalar(0.0);
    cmsOnCmsMatrix_->PutScalar(0.0);
    densityOnDensityMatrix_->PutScalar(0.0);
    densityOnCmsMatrix_->PutScalar(0.0);
    poissonMatrix_->PutScalar(0.0);
    cmsOnPoissonMatrix_->PutScalar(0.0);
    poissonOnDensityMatrix_->PutScalar(0.0);
  }
  return(0);
}
//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {

  if (poissonMap_.MyGID(rowGID)) {
    if (poissonMap_.MyGID(colGID)) {
      if (firstTime_) {
	if (rowGID != curPoissonRow_) {
	  insertPoissonRow();
	  curPoissonRow_ = rowGID;
	}
	curPoissonRowValues_[colGID] += value;
      }
      else {
	poissonMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
      }
    }
    else { //!poissonMap_.MyGID(colGID)
      if (firstTime_) {
	if (rowGID != curPDRow_) {
	  insertPDRow();
	  curPDRow_ = rowGID;
	}
	curPDRowValues_[colGID] += value;
      }
      else {
	poissonOnDensityMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
      }
    }
  } //end if poissonMap_.MyGID(rowGID)

  else if (cmsMap_.MyGID(rowGID)) {
    if (rowGID==colGID)
      (*cmsOnCmsMatrix_)[cmsMap_.LID(rowGID)] += value; // Storing this cms block in a vector since it is diagonal
    else if (densityMap_.MyGID(colGID)){
      if (firstTime_) {
	if (rowGID!=curRow_) { 
	  insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
	  curRow_=rowGID;
	}
	curRowValues_[colGID] += value;
      }
      else if (!isFLinear_) {
	int newRowGID = densityMap_.GID(cmsMap_.LID(rowGID));
	cmsOnDensityMatrix_->SumIntoGlobalValues(newRowGID, 1, &value, &colGID);
      }
    }
    else if (poissonMap_.MyGID(colGID)) {
      if (firstTime_) {
	if (rowGID!=curCPRow_) {
	  insertCPRow(); // Dump the current contents of curCPRowValues_ into matrix and clear map
	  curCPRow_=rowGID;
	}
	curCPRowValues_[colGID] += value;
      }
      else {
	cmsOnPoissonMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
      }
    }
  } //end if cmsMap_.MyGID(rowGID)

  else {
    if (rowGID==colGID)
      (*densityOnDensityMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
    else
      (*densityOnCmsMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
  }

  return(0);
}
//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::insertPoissonRow() {

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
int dft_PolyA22_Coulomb_Epetra_Operator::insertCPRow() {

  if (curCPRowValues_.empty()) return(0);
  int numEntries = curCPRowValues_.size();
  if (numEntries>indices_.Length()) {
    indices_.Resize(numEntries);
    values_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos=curCPRowValues_.begin(); pos!=curCPRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  cmsOnPoissonMatrix_->InsertGlobalValues(curCPRow_, numEntries, values_.Values(), indices_.Values());
  curCPRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::insertPDRow() {

  if (curPDRowValues_.empty()) return(0);
  int numEntries = curPDRowValues_.size();
  if (numEntries>indices_.Length()) {
    indices_.Resize(numEntries);
    values_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos=curPDRowValues_.begin(); pos!=curPDRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  poissonOnDensityMatrix_->InsertGlobalValues(curPDRow_, numEntries, values_.Values(), indices_.Values());
  curPDRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::finalizeProblemValues() {

  if (isLinearProblemSet_) return(0); // nothing to do

  if (!isFLinear_) {
    insertRow(); // Dump any remaining entries
    cmsOnDensityMatrix_->FillComplete();
    cmsOnDensityMatrix_->OptimizeStorage();
  }

  if (firstTime_) {
    insertPoissonRow();
    insertCPRow();
    insertPDRow();
  }

  poissonMatrix_->FillComplete();
  poissonMatrix_->OptimizeStorage();

  cmsOnPoissonMatrix_->FillComplete(poissonMap_, cmsMap_);
  cmsOnPoissonMatrix_->OptimizeStorage();

  poissonOnDensityMatrix_->FillComplete(densityMap_, poissonMap_);
  poissonOnDensityMatrix_->OptimizeStorage();

  isLinearProblemSet_ = true;
  firstTime_ = false;

  return(0);
}
//==============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  // The true A22 block is of the form:

  // |  Dcc     F    |
  // |  Ddc     Ddd  |
  
  // where 
  // Dcc is Cms on Cms (diagonal),
  // F is Cms on Density (fairly dense)
  // Ddc is Density on Cms (diagonal with small coefficient values),
  // Ddd is Density on Density (diagonal).
  //
  // We will approximate A22 with:
  
  // |  Dcc     F    |
  // |  0       Ddd  |
  
  // replacing Ddc with a zero matrix for the ApplyInverse method only.

  // Our algorithm is then:
  // Y2 = Ddd \ X2
  // Y1 = Dcc \ (X1 - F*Y2)

  // If F is in SW:
  // The true A22 block is of the form:

  // |  Ddd     Ddc  |
  // |  F       Dcc  |
  
  // where 
  // Ddd is Density on Density (diagonal),
  // Ddc is Density on Cms (diagonal with small coefficient values),
  // F is Cms on Density (fairly dense),
  // Dcc is Cms on Cms (diagonal).
  //
  // We will approximate A22 with:
  
  // |  Ddd     0    |
  // |  F       Dcc  |
  
  // replacing Ddc with a zero matrix for the ApplyInverse method only.

  // Our algorithm is then:
  // Y1 = Ddd \ X1
  // Y2 = Dcc \ (X2 - F*Y1)  

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numCmsElements = cmsMap_.NumMyElements();
  int numDensityElements = densityMap_.NumMyElements();
  int numPoissonElements = poissonMap_.NumMyElements();

  double ** X0ptr;
  double ** Y0ptr;
  double ** X1ptr = new double *[NumVectors];
  double ** Y1ptr = new double *[NumVectors];
  double ** X2ptr = new double *[NumVectors];
  double ** Y2ptr = new double *[NumVectors];

  X.ExtractView(&X0ptr); // Get array of pointers to columns of X
    Y.ExtractView(&Y0ptr); // Get array of pointers to columns of Y
  if (F_location_ == 1) { //F in NE part
    for (int i=0; i<NumVectors; i++) {
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numCmsElements;
      X2ptr[i] = X1ptr[i]+numCmsElements;
    }
  }
  else { //F in SW part
    for (int i = 0;i<NumVectors; i++) {
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      X2ptr[i] = X1ptr[i]+numDensityElements;
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numDensityElements;
    }
  }
  Epetra_MultiVector X0(View, poissonMap_, X0ptr, NumVectors); // X0 is a view of the first numPoisson elements of X
  Epetra_MultiVector X1(View, densityMap_, X1ptr, NumVectors); // X1 is a view of the middle numDensity/numCms elements of X
  Epetra_MultiVector X2(View, densityMap_, X2ptr, NumVectors); // X2 is a view of the last numDensity/numCms elements of X
  Epetra_MultiVector Y0(View, poissonMap_, Y0ptr, NumVectors); // Y0 is a view of the first numPoisson elements of Y
  Epetra_MultiVector Y1(View, densityMap_, Y1ptr, NumVectors); // Y1 is a view of the middle numDensity/numCms elements of Y
  Epetra_MultiVector Y2(View, densityMap_, Y2ptr, NumVectors); // Y2 is a view of the last numDensity/numCms elements of Y

  Epetra_MultiVector Y0tmp(Y0);
  Epetra_MultiVector Y1tmp(Y1);
  Epetra_MultiVector Y1tmp2(Y1);
  Epetra_MultiVector Y2tmp(Y2);
  Epetra_MultiVector Y2tmp2(Y2);

  if (F_location_ == 1) {
    Y2.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X2, 0.0);
    poissonOnDensityMatrix_->Apply(Y2, Y0tmp);
    Y0.Update(1.0, X0, -1.0, Y0tmp, 0.0);
  }
  else {
    Y1.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X1, 0.0);
    poissonOnDensityMatrix_->Apply(Y1, Y0tmp);
    Y0.Update(1.0, X0, -1.0, Y0tmp, 0.0);
  }

  //now apply inverse of poissonMatrix_ to Y0
  Epetra_LinearProblem implicitProblem(&(*poissonMatrix_), &Y0, &Y0);
  Teuchos::RefCountPtr<AztecOO> solver = Teuchos::rcp(new AztecOO(implicitProblem));
  ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(*poissonMatrix_, MLList_);
  
  solver->SetParameters(*parameterList_);
  //if (options_!=0) solver->SetAllAztecOptions(options_);
  // if (params_!=0) solver->SetAllAztecParams(params_);
  
  //const int * options = solver->GetAllAztecOptions();
  //const double * params = solver->GetAllAztecParams();
  solver->SetPrecOperator(MLPrec);
  //solver->SetAztecOption(AZ_solver, AZ_cg); //change?
  solver->SetAztecOption(AZ_solver, AZ_gmres);
  solver->SetAztecOption(AZ_output, 1); //change?
  solver->Iterate(Teuchos::getParameter<int>(*parameterList_, "Max_iter"), Teuchos::getParameter<double>(*parameterList_, "Tol")); //Try to solve
  //solver->Iterate(options[AZ_max_iter], params[AZ_tol]); // Try to solve
  delete MLPrec;
  
  if (F_location_ == 1) {
    cmsOnDensityMatrix_->Apply(Y2, Y1tmp);
    cmsOnPoissonMatrix_->Apply(Y0, Y1tmp2);
    Y1tmp.Update(1.0, Y1tmp2, 1.0);
    Y1.Update(1.0, X1, -1.0, Y1tmp, 0.0);
    Y1.ReciprocalMultiply(1.0, *cmsOnCmsMatrix_, Y1, 0.0);
  }
  else { 
    cmsOnDensityMatrix_->Apply(Y1, Y2tmp);
    cmsOnPoissonMatrix_->Apply(Y0, Y2tmp2);
    Y2tmp.Update(1.0, Y2tmp2, 1.0);
    Y2.Update(1.0, X2, -1.0, Y2tmp, 0.0);
    Y2.ReciprocalMultiply(1.0, *cmsOnCmsMatrix_, Y2, 0.0);
  }
  
  delete [] Y1ptr;
  delete [] X1ptr;
  delete [] Y2ptr;
  delete [] X2ptr;
    
  return(0);
}
//==============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numCmsElements = cmsMap_.NumMyElements();
  int numDensityElements = densityMap_.NumMyElements();
  int numPoissonElements = poissonMap_.NumMyElements();

  double ** X0ptr;
  double ** Y0ptr;
  double ** X1ptr = new double *[NumVectors];
  double ** Y1ptr = new double *[NumVectors];
  double ** X2ptr = new double *[NumVectors];
  double ** Y2ptr = new double *[NumVectors];

  X.ExtractView(&X0ptr); // Get array of pointers to columns of X
  Y.ExtractView(&Y0ptr); // Get array of pointers to columns of Y
  if (F_location_ == 1) {
    for (int i=0; i<NumVectors; i++) {
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      X2ptr[i] = X1ptr[i]+numCmsElements;
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numCmsElements;
    }
  }
  else {
    for (int i=0; i<NumVectors; i++) {
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      X2ptr[i] = X1ptr[i]+numDensityElements;
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numDensityElements;
    }
  }

 Epetra_MultiVector X0(View, poissonMap_, X0ptr, NumVectors); // Start X0 to view the first numPoisson elements of X
  Epetra_MultiVector X1(View, densityMap_, X1ptr, NumVectors); // Start X1 to view middle numDensityElements/numCms elements of X
  Epetra_MultiVector X2(View, densityMap_, X2ptr, NumVectors); // Start X2 to view last numDensity/numCms elements of X - was cms
  Epetra_MultiVector Y0(View, poissonMap_, Y0ptr, NumVectors); // Y0 is a view of the first numPoisson elements of Y
  Epetra_MultiVector Y1(View, densityMap_, Y1ptr, NumVectors); // Y1 is a view of the middle numDensity/numCms elements of Y
  Epetra_MultiVector Y2(View, densityMap_, Y2ptr, NumVectors); // Y2 is a view of the last numDensity/numCms elements of Y
  Epetra_MultiVector Y0tmp(Y0);
  Epetra_MultiVector Y0tmp2(Y0);
  Epetra_MultiVector Y1tmp(Y1);
  Epetra_MultiVector Y1tmp2(Y1);
  Epetra_MultiVector Y2tmp(Y2);
  Epetra_MultiVector Y2tmp2(Y2);

  poissonMatrix_->Apply(X0, Y0tmp);

  if (F_location_ == 1) {
    poissonOnDensityMatrix_->Apply(X2, Y0tmp2);
    cmsOnPoissonMatrix_->Apply(X0, Y1tmp);
    cmsOnDensityMatrix_->Apply(X2, Y1tmp2);
    Y2tmp.Multiply(1.0, *densityOnCmsMatrix_, X1, 0.0);
    Y2tmp2.Multiply(1.0, *densityOnDensityMatrix_, X2, 0.0);
    Y1.Multiply(1.0, *cmsOnCmsMatrix_, X1, 0.0);
    Y1.Update(1.0, Y1tmp, 1.0,  Y1tmp2, 1.0);
    Y2.Update(1.0, Y2tmp, 1.0, Y2tmp2, 0.0);
    Y0.Update(1.0, Y0tmp, 1.0, Y0tmp2, 0.0);
  }
  else {
    poissonOnDensityMatrix_->Apply(X1, Y0tmp2);
    Y1tmp.Multiply(1.0, *densityOnDensityMatrix_, X1, 0.0);
    Y1tmp2.Multiply(1.0, *densityOnCmsMatrix_, X2, 0.0);
    cmsOnPoissonMatrix_->Apply(X0, Y2tmp);
    cmsOnDensityMatrix_->Apply(X1, Y2tmp2);
    Y2.Multiply(1.0, *cmsOnCmsMatrix_, X2, 0.0);
    Y2.Update(1.0, Y2tmp, 1.0, Y2tmp2, 1.0);
    Y1.Update(1.0, Y1tmp, 1.0, Y1tmp2, 0.0);
    Y0.Update(1.0, Y0tmp, 1.0, Y0tmp2, 0.0);
  }
  
  delete [] X1ptr;
  delete [] Y1ptr;
  delete [] X2ptr;
  delete [] Y2ptr;
  
  return(0);
}
//==============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::Check(bool verbose) const {

  Epetra_Vector x(OperatorDomainMap());
  Epetra_Vector b(OperatorRangeMap());
  x.Random(); // Fill x with random numbers

  Epetra_Vector x1(View, densityMap_, x.Values()+poissonMap_.NumMyElements()); // Start x1 to view middle numCmsElements elements of x
  Epetra_Vector b2(View, densityMap_, b.Values()+poissonMap_.NumMyElements()+cmsMap_.NumMyElements()); // Start b2 to view last numDensity elements of b

  //SMK
  Epetra_Vector x2(View, densityMap_, x.Values()+poissonMap_.NumMyElements()+densityMap_.NumMyElements()); //Start x2 to view last numCms elements of x
  Epetra_Vector b1(View, densityMap_, b.Values()+poissonMap_.NumMyElements()); // Start b1 to view middle numDensity elements of b
  
  Apply(x, b); // Forward operation
  
  if (F_location_ == 1) {
    // Inverse is not exact, so we must modify b2 first:
    b2.Multiply(-1.0, *densityOnCmsMatrix_, x1, 1.0);
  }
  else {
    // Inverse is not exact, so we must modify b1 first:
    b1.Multiply(-1.0, *densityOnCmsMatrix_, x2, 1.0);
  }
  
  ApplyInverse(b, b); // Reverse operation
  b.Update(-1.0, x, 1.0); // Should be zero
  double resid = 0.0;
  b.Norm2(&resid);
  
  if (verbose) 
    std::cout << "A22 self-check residual = " << resid << endl;

  if (resid > 1.0E-12) return(-1); // Bad residual
  return(0);
}
