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

#include "sol_A22_G11_Z12_D21_D22.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Comm.h"
#include "Epetra_Distributor.h"
#include "EpetraExt_RowMatrixOut.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Teuchos_TestForException.hpp"

//=============================================================================
sol_A22_G11_Z12_D21_D22::sol_A22_G11_Z12_D21_D22(const Epetra_Map & map1, const Epetra_Map & map2, 
                                                 const Epetra_Map & A22Map, Teuchos::ParameterList * parameterList) 
  : map1_(map1),
    map2_(map2),
    A22Map_(A22Map),
    parameterList_(parameterList),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1) {

      G11Matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, map1, 0));
      D21Matrix_ = Teuchos::rcp(new Epetra_Vector(map2));
      D22Matrix_ = Teuchos::rcp(new Epetra_Vector(map2));
      Label_ = "sol_A22_G11_Z12_D21_D22 operator";
      G11Matrix_->SetLabel("A22_G11_Matrix");
      isG11RowStatic_.resize(map1.NumMyElements());
      isG11RowStatic_.assign(map1.NumMyElements(), false); //  all rows initialized as non-static

  }
//==============================================================================
sol_A22_G11_Z12_D21_D22::~sol_A22_G11_Z12_D21_D22() {
}
//=============================================================================
int sol_A22_G11_Z12_D21_D22::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    D21Matrix_->PutScalar(0.0);
    D22Matrix_->PutScalar(0.0);
    int numEntries;
    double * values;
    for (int i = 0 ; i < G11Matrix_->NumMyRows() ; ++i) {
      if (!(isG11RowStatic_[i])) {
        G11Matrix_->ExtractMyRowView(i, numEntries, values);
        for (int k = 0; k < numEntries; ++k) values[k] = 0.0;
      }
    }
    
  }
  
  return(0);
}
//=============================================================================
int sol_A22_G11_Z12_D21_D22::insertMatrixValue(int rowGID, int colGID, double value) {


  if (map1_.MyGID(rowGID)) {   // G11
    TEST_FOR_EXCEPTION(map2_.MyGID(colGID), std::logic_error, 
                       "You are adding a value to the 12 block of A22.  This action is not allowed since the 12 block by definition zero"); 
    if (firstTime_) {
      if (rowGID!=curRow_) { 
        insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
        curRow_=rowGID;
      }
      curRowValues_[colGID] += value;
    }
    else { 
      TEST_FOR_EXCEPTION(isG11RowStatic_[map1_.LID(rowGID)], std::logic_error, 
                         "You are adding a value to the 21 block of A22.  This action is not allowed since the 21 block by definition zero"); 
      G11Matrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
    }
  }
  else { // D21, D22
    TEST_FOR_EXCEPTION(!map2_.MyGID(rowGID), std::logic_error, 
                       "You are adding values to a row that is not part of A22"); 
    if (rowGID==colGID)
      (*D22Matrix_)[map2_.LID(rowGID)] += value; // Storing this block in a vector since it is diagonal
    else {
      TEST_FOR_EXCEPT(map2_.LID(rowGID)!=map1_.LID(colGID)); // Confirm that this is a diagonal value
      (*D21Matrix_)[map2_.LID(rowGID)] += value; // Storing this block in a vector since it is diagonal
    }
  }
  
  return(0);
}
//=============================================================================
int sol_A22_G11_Z12_D21_D22::insertRow() {
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
  
  G11Matrix_->InsertGlobalValues(curRow_, numEntries, values_.Values(), indices_.Values());

  curRowValues_.clear();
  return(0);
}
//=============================================================================
int sol_A22_G11_Z12_D21_D22::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  insertRow(); // Dump any remaining entries
  G11Matrix_->FillComplete(map1_, map1_);

  if (firstTime_) {
    Teuchos::ParameterList list;
    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    string PrecType = "ILU"; // incomplete LU
    //string PrecType = "point relaxation"; // Gauss-Seidel
    int OverlapLevel = 0; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.
    
    G11Inverse_  = Teuchos::rcp(factory_.Create(PrecType, G11Matrix_.get(), OverlapLevel));
    assert(G11Inverse_.get()!=0);
    
    // specify parameters for ILU
    list.set("fact: drop tolerance", 1e-9);  // these should be input parameters from Tramonto
    list.set("fact: level-of-fill", 2);
    //list.set("relaxation: type", "Gauss-Seidel");
    
    // sets the parameters
    IFPACK_CHK_ERR(G11Inverse_->SetParameters(list));
    
    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    IFPACK_CHK_ERR(G11Inverse_->Initialize());
  }
  // Builds the preconditioners, by looking for the values of
  // the matrix.
  IFPACK_CHK_ERR(G11Inverse_->Compute());
  
  isLinearProblemSet_ = true;
  firstTime_ = false;

  return(0);
}
//==============================================================================
int sol_A22_G11_Z12_D21_D22::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  
  // Solve Y = inv(A22)*X where X is known and inv(A22) is an approximation to the inverse of A22.
  
  // The  A22 block is of the form:
  
  // |  Gff      0   |
  // |  Dfd      Ddd |
  
  // where 
  // Gff is field on field (general),
  // Dfd is field on Density 
  // Ddd is Density on Density (diagonal).
  //
  
  // The equations are 
  // Gff*Y1          = X1
  // Dfd*Y1 + Ddd*Y2 = X2
  
  // Our algorithm is then:
  // Y1 =Ifpack(Gff) \ X1
  // Y2 = Ddd \ (X2 - Dfd*Y1)
  
  // where Ifpack(Gff) is a preconditioner selected via the solver parameters using Ifpack methods.  
  
  
  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numElements1 = map1_.NumMyElements();
  int numElements2 = map2_.NumMyElements();
  
  double ** Y1ptr;
  double ** X1ptr;
  double ** Y2ptr = new double *[NumVectors];
  double ** X2ptr = new double *[NumVectors];
  
  Y.ExtractView(&Y1ptr); // Get array of pointers to columns of Y
  X.ExtractView(&X1ptr); // Get array of pointers to columns of X
  for (int i=0; i<NumVectors; i++) {
    Y2ptr[i] = Y1ptr[i]+numElements1;
    X2ptr[i] = X1ptr[i]+numElements1;
  }
  Epetra_MultiVector Y1(View, map1_, Y1ptr, NumVectors); // Y1 is a view of the first numField elements of Y
  Epetra_MultiVector Y2(View, map2_, Y2ptr, NumVectors); // Start Y2 to view last numDensity elements of Y
  Epetra_MultiVector X1(View, map1_, X1ptr, NumVectors); // Start X1 to view first numField elements of X
  Epetra_MultiVector X2(View, map2_, X2ptr, NumVectors); // Start X2 to view last numDensity elements of X
  
  Epetra_MultiVector X2tmp(X2);
  TEST_FOR_EXCEPT(G11Inverse_->ApplyInverse(X1,Y1)!=0);                      // Ifpack solve Y1 = G11 \ X1
  TEST_FOR_EXCEPT(X2tmp.Multiply(-11.0, *D21Matrix_, Y1, 1.0)!=0);           // X2tmp = (X2 - Dfd*Y1)
  TEST_FOR_EXCEPT(Y2.ReciprocalMultiply(1.0, *D22Matrix_, X2tmp, 0.0)!=0);   // Y2 = Ddd \ X2tmp
  
  delete [] Y2ptr;
  delete [] X2ptr;
  
  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA22 ApplyInverse = " << normvalue << endl;
  return(0);
}
//==============================================================================
int sol_A22_G11_Z12_D21_D22::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  //double normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in PolyA22 Apply = " << normvalue << endl;

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numElements1 = map1_.NumMyElements();
  int numElements2 = map2_.NumMyElements(); 

  // D21Matrix will be nonzero only if field and density maps are the same size
  bool hasDensityOnField = map1_.NumGlobalElements()==map2_.NumGlobalElements(); 

  double ** X1ptr;
  double ** Y1ptr;
  double ** X2ptr = new double *[NumVectors];
  double ** Y2ptr = new double *[NumVectors];

  X.ExtractView(&X1ptr); // Get array of pointers to columns of X
  Y.ExtractView(&Y1ptr); // Get array of pointers to columns of Y
  for (int i=0; i<NumVectors; i++) {
    X2ptr[i] = X1ptr[i]+numElements1;
    Y2ptr[i] = Y1ptr[i]+numElements1;
  }
  
  Epetra_MultiVector Y1(View, map1_, Y1ptr, NumVectors); // Y1 is a view of the first numField elements of Y
  Epetra_MultiVector Y2(View, map2_, Y2ptr, NumVectors); // Start Y2 to view last numDensity elements of Y
  Epetra_MultiVector X1(View, map1_, X1ptr, NumVectors); // Start X1 to view first numField elements of X
  Epetra_MultiVector X2(View, map2_, X2ptr, NumVectors); // Start X2 to view last numDensity elements of X
  
  G11Matrix_->Apply(X1, Y1);
  TEST_FOR_EXCEPT(Y2.Multiply(1.0, *D21Matrix_, X1, 0.0)!=0);
  TEST_FOR_EXCEPT(Y2.Multiply(1.0, *D22Matrix_, X2, 1.0)!=0);
  
  delete [] X2ptr;
  delete [] Y2ptr;

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA22 Apply = " << normvalue << endl;
  return(0);
}
//==============================================================================
int sol_A22_G11_Z12_D21_D22::Check(bool verbose) const {

  Epetra_Vector x(OperatorDomainMap());
  Epetra_Vector b(OperatorRangeMap());
  x.Random(); // Fill x with random numbers

  Apply(x, b); // Forward operation

  // Inverse is not exact, so we must modify b1 first:
  // b1 = b1 - (G11*x1 - Ifpack(G11)*x1)
  Epetra_Vector x1(View, map1_, x.Values()); // Start x1, b1 to view first numElements1 elements of x
  Epetra_Vector G11x1(map1_);
  Epetra_Vector IfpackG11x1(map1_);
  TEST_FOR_EXCEPT(G11Matrix_->Apply(x1,G11x1)!=0);
  TEST_FOR_EXCEPT(G11Inverse_->Apply(x1,IfpackG11x1)!=0);
  Epetra_Vector b1(View, map1_, b.Values()); 
  b1.Update(-1.0, G11x1, 1.0, IfpackG11x1, 1.0);

  ApplyInverse(b, b); // Reverse operation

  b.Update(-1.0, x, 1.0); // Should be zero

  double resid = 0.0;
  b.Norm2(&resid);

  if (verbose) 
    std::cout << "A22 self-check residual = " << resid << endl;

  if (resid > 1.0E-12) return(-1); // Bad residual
  return(0);
}
