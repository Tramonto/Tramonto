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

#include "dft_BasicA22_Epetra_Operator.hpp"
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
dft_BasicA22_Epetra_Operator::dft_BasicA22_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & block2Map, Teuchos::ParameterList * parameterList) 
  : cmsMap_(cmsMap),
    densityMap_(densityMap),
    block2Map_(block2Map),
    parameterList_(parameterList),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    isFLinear_(false),
    firstTime_(true),
    curRow_(-1) {

  A22Matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, block2Map, 0));
  isStaticRow_.resize(block2Map.NumMyElements()); 
  isStaticRow_.assign(block2Map.NumMyElements(), false); //  all rows initialized as non-static
  Label_ = "dft_BasicA22_Epetra_Operator";
  cmsOnDensityMatrix_->SetLabel("BasicA22::A22Matrix");
  F_location_ = Teuchos::getParameter<int>(*parameterList_, "F_location"); //F in NE if F_location = 1, F in SW otherwise
  }
//==============================================================================
dft_BasicA22_Epetra_Operator::~dft_BasicA22_Epetra_Operator() {
}
//=============================================================================
int dft_BasicA22_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    int numEntries;
    double * values;
    for (int i = 0 ; i < A22Matrix_->NumMyRows() ; ++i) {
      if (!(isStaticRow_[i])) {
        A22Matrix_->ExtractMyRowView(i, numEntries, values);
        for (int k = 0; k < numEntries; ++k) values[k] = 0.0;
		  }
    }
  }
  
  return(0);
}
//=============================================================================
int dft_BasicA22_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {

  
  if (cmsMap_.MyGID(rowGID)) {
    if (rowGID==colGID)
      (*cmsOnCmsMatrix_)[cmsMap_.LID(rowGID)] += value; // Storing this cms block in a vector since it is diagonal
    else {
      if (firstTime_) {
	if (rowGID!=curRow_) { 
	  insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
	  curRow_=rowGID;
	}
	curRowValues_[colGID] += value;
      }
      else if (!isFLinear_) { 
	cmsOnDensityMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
      }
    }
  }
  else {
    if (rowGID==colGID)
      (*densityOnDensityMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
    else {
      TEST_FOR_EXCEPT(densityMap_.LID(rowGID)!=cmsMap_.LID(colGID)); // Confirm that this is a diagonal value
      (*densityOnCmsMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
    }
  }

  return(0);
}
//=============================================================================
int dft_BasicA22_Epetra_Operator::insertRow() {
  if (curRowValues_.empty()) return(0);
  int numEntries = curRowValues_.size();
  if (numEntries>indices_.Length()) {
    indices_.Resize(numEntries);
    values_.Resize(numEntries);
  }
  int i=0;
  std::map<int, double>::iterator pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos) {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  }
  
  cmsOnDensityMatrix_->InsertGlobalValues(curRow_, numEntries, values_.Values(), indices_.Values());

  curRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_BasicA22_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do
  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_.NumGlobalElements()==densityMap_.NumGlobalElements(); 

  if (!isFLinear_) {
    insertRow(); // Dump any remaining entries
    cmsOnDensityMatrix_->FillComplete(densityMap_, cmsMap_);
    cmsOnDensityMatrix_->OptimizeStorage();
  }

  if (!hasDensityOnCms) { // Confirm that densityOnCmsMatrix is zero
    double normvalue;
    densityOnCmsMatrix_->NormInf(&normvalue);
    TEST_FOR_EXCEPT(normvalue!=0.0);
  }
  //cout << "CmsOnDensityMatrix Inf Norm = " << cmsOnDensityMatrix_->NormInf() << endl;
  //densityOnDensityMatrix_->NormInf(&normvalue);
  //cout << "DensityOnDensityMatrix Inf Norm = " << normvalue << endl;
  //cout << "DensityOnCmsMatrix Inf Norm = " << normvalue << endl;

  isLinearProblemSet_ = true;
  firstTime_ = false;

  return(0);
}
//==============================================================================
int dft_BasicA22_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
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

  // Or, if F is in the SW quadrant:
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

  //double normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in BasicA22 ApplyInverse = " << normvalue << endl;

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numCmsElements = cmsMap_.NumMyElements();
  int numDensityElements = densityMap_.NumMyElements();

  double ** Y1ptr;
  double ** X1ptr;
  double ** Y2ptr = new double *[NumVectors];
  double ** X2ptr = new double *[NumVectors];

  Y.ExtractView(&Y1ptr); // Get array of pointers to columns of Y
  X.ExtractView(&X1ptr); // Get array of pointers to columns of X
  if (F_location_ == 1) { //F in NE
    for (int i=0; i<NumVectors; i++) {
      Y2ptr[i] = Y1ptr[i]+numCmsElements;
      X2ptr[i] = X1ptr[i]+numCmsElements;
    }
    Epetra_MultiVector Y1(View, cmsMap_, Y1ptr, NumVectors); // Y1 is a view of the first numCms elements of Y
    Epetra_MultiVector Y2(View, densityMap_, Y2ptr, NumVectors); // Start Y2 to view last numDensity elements of Y
    Epetra_MultiVector X1(View, cmsMap_, X1ptr, NumVectors); // Start X1 to view first numCms elements of X
    Epetra_MultiVector X2(View, densityMap_, X2ptr, NumVectors); // Start X2 to view last numDensity elements of X
    
    Epetra_MultiVector Y1tmp(Y1);
    TEST_FOR_EXCEPT(Y2.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X2, 0.0)!=0);
    TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(Y2, Y1tmp)!=0);
    TEST_FOR_EXCEPT(Y1.Update(1.0, X1, -1.0, Y1tmp, 0.0)!=0);
    TEST_FOR_EXCEPT(Y1.ReciprocalMultiply(1.0, *cmsOnCmsMatrix_, Y1, 0.0)!=0);
  }
  else {
    for (int i = 0;i<NumVectors; i++) {
      Y2ptr[i] = Y1ptr[i]+numDensityElements;
      X2ptr[i] = X1ptr[i]+numDensityElements;
    }
    
    Epetra_MultiVector Y1(View, densityMap_, Y1ptr, NumVectors); // Y1 is a view of the first numDensity elements of Y
    Epetra_MultiVector Y2(View, cmsMap_, Y2ptr, NumVectors); // Start Y2 to view last numCms elements of Y
    Epetra_MultiVector X1(View, densityMap_, X1ptr, NumVectors); // Start X1 to view first numDensity elements of X
    Epetra_MultiVector X2(View, cmsMap_, X2ptr, NumVectors); // Start X2 to view last numCms elements of X
    Epetra_MultiVector Y2tmp(Y2);
    TEST_FOR_EXCEPT(Y1.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X1, 0.0)!=0);
    TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(Y1, Y2tmp)!=0);
    TEST_FOR_EXCEPT(Y2.Update(1.0, X2, -1.0, Y2tmp, 0.0)!=0);
    TEST_FOR_EXCEPT(Y2.ReciprocalMultiply(1.0, *cmsOnCmsMatrix_, Y2, 0.0)!=0);
  }
  
  delete [] Y2ptr;
  delete [] X2ptr;

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in BasicA22 ApplyInverse = " << normvalue << endl;
  return(0);
}
//==============================================================================
int dft_BasicA22_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  //double normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in BasicA22 Apply = " << normvalue << endl;

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numCmsElements = cmsMap_.NumMyElements();
  int numDensityElements = densityMap_.NumMyElements(); 

  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_.NumGlobalElements()==densityMap_.NumGlobalElements(); 

  double ** X1ptr;
  double ** Y1ptr;
  double ** X2ptr = new double *[NumVectors];
  double ** Y2ptr = new double *[NumVectors];

  X.ExtractView(&X1ptr); // Get array of pointers to columns of X
  Y.ExtractView(&Y1ptr); // Get array of pointers to columns of Y
  if (F_location_ == 1) {//F in NE
    for (int i=0; i<NumVectors; i++) {
      X2ptr[i] = X1ptr[i]+numCmsElements;
      Y2ptr[i] = Y1ptr[i]+numCmsElements;
    }
  
    Epetra_MultiVector Y1(View, cmsMap_, Y1ptr, NumVectors); // Y1 is a view of the first numCms elements of Y
    Epetra_MultiVector Y2(View, densityMap_, Y2ptr, NumVectors); // Start Y2 to view last numDensity elements of Y
    Epetra_MultiVector X1(View, cmsMap_, X1ptr, NumVectors); // Start X1 to view first numCms elements of X
    Epetra_MultiVector X2(View, densityMap_, X2ptr, NumVectors); // Start X2 to view last numDensity elements of X

    cmsOnDensityMatrix_->Apply(X2, Y1);
    TEST_FOR_EXCEPT(Y1.Multiply(1.0, *cmsOnCmsMatrix_, X1, 1.0)!=0);
    TEST_FOR_EXCEPT(Y2.Multiply(1.0, *densityOnCmsMatrix_, X1, 0.0)!=0);
    TEST_FOR_EXCEPT(Y2.Multiply(1.0, *densityOnDensityMatrix_, X2, 1.0)!=0);
  }
  else {
    for (int i=0; i<NumVectors; i++) {
      X2ptr[i] = X1ptr[i]+numDensityElements;
      Y2ptr[i] = Y1ptr[i]+numDensityElements;
    }
  
    Epetra_MultiVector Y1(View, densityMap_, Y1ptr, NumVectors); // Y1 is a view of the first numDensity elements of Y
    Epetra_MultiVector Y2(View, cmsMap_, Y2ptr, NumVectors); // Start Y2 to view last numCms elements of Y
    Epetra_MultiVector X1(View, densityMap_, X1ptr, NumVectors); // Start X1 to view first numDensity elements of X
    Epetra_MultiVector X2(View, cmsMap_, X2ptr, NumVectors); // Start X2 to view last numCms elements of X
    if (hasDensityOnCms) {
      // convert X2 map
      Epetra_MultiVector X2tmp(View, densityMap_, X2ptr, NumVectors);
      TEST_FOR_EXCEPT(Y1.Multiply(1.0, *densityOnCmsMatrix_, X2tmp, 0.0)!=0); // Momentarily make X2 compatible with densityMap_
    }
    TEST_FOR_EXCEPT(Y1.Multiply(1.0, *densityOnDensityMatrix_, X1, 1.0)!=0);
    TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(X1, Y2)!=0);
    TEST_FOR_EXCEPT(Y2.Multiply(1.0, *cmsOnCmsMatrix_, X2, 1.0)!=0);
  }
  
  delete [] X2ptr;
  delete [] Y2ptr;

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in BasicA22 Apply = " << normvalue << endl;
  return(0);
}
//==============================================================================
int dft_BasicA22_Epetra_Operator::Check(bool verbose) const {

  Epetra_Vector x(OperatorDomainMap());
  Epetra_Vector b(OperatorRangeMap());
  x.Random(); // Fill x with random numbers

  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_.NumGlobalElements()==densityMap_.NumGlobalElements(); 

  Apply(x, b); // Forward operation

  if (hasDensityOnCms) {
      
    if (F_location_ == 1) {//F in NE
      // Inverse is not exact, so we must modify b2 first:
      Epetra_Vector x1(View, densityMap_, x.Values()); // Start x1 to view first numCmsElements elements of x
      Epetra_Vector b2(View, densityMap_, b.Values()+cmsMap_.NumMyElements()); // Start b2 to view last numDensity elements of b
      b2.Multiply(-1.0, *densityOnCmsMatrix_, x1, 1.0);
    }
    else {
      // Inverse is not exact, so we must modify b1 first:
      Epetra_Vector x2(View, densityMap_, x.Values()+densityMap_.NumMyElements()); //Start x2 to view last numCms elements of x
      Epetra_Vector b1(View, densityMap_, b.Values()); // Start b1 to view first numDensity elements of b
      b1.Multiply(-1.0, *densityOnCmsMatrix_, x2, 1.0);
    }
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
