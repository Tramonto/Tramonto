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

#include "dft_Poisson_Epetra_Operator.hpp"
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

//==============================================================================
dft_Poisson_Epetra_Operator::dft_Poisson_Epetra_Operator(const Epetra_Map & poissonMap) 
  : poissonMap_(poissonMap),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true),
    curRow_(-1) {

  matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, poissonMap, 0));
  Label_ = "dft_Poisson_Epetra_Operator";
  matrix_->SetLabel("dft_Poisson_Epetra_Operator::Poisson");
}
//==============================================================================
dft_Poisson_Epetra_Operator::~dft_Poisson_Epetra_Operator() {
}
//=============================================================================
int dft_Poisson_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    matrix_->PutScalar(0.0);
  
  return(0);
}
//=============================================================================
int dft_Poisson_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {
  
  
  if (firstTime_) {
    if (rowGID!=curRow_) { 
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    }
    curRowValues_[colGID] += value;
  }
  else
    matrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);

  return(0);
}
//=============================================================================
int dft_Poisson_Epetra_Operator::insertRow() {
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
  matrix_->InsertGlobalValues(curRow_, numEntries, values_.Values(), indices_.Values());

  curRowValues_.clear();
  return(0);
}
//=============================================================================
int dft_Poisson_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  insertRow(); // Dump any remaining entries
  matrix_->FillComplete();
  matrix_->OptimizeStorage();
  //EpetraExt::RowMatrixToMatrixMarketFile( "Poisson.dat", matrix_, "Poisson operator blocks", 
  //				  "The Poisson block of  problem",
  //				  true);
 // abort();
  
  /*  std::cout << matrix_<< std::endl;
  */

  if (firstTime_) {
    Teuchos::ParameterList list;
    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    string PrecType = "ILU"; // incomplete LU
    //string PrecType = "point relaxation"; // Gauss-Seidel
    int OverlapLevel = 0; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.
    
    ifpackInverse_  = Teuchos::rcp(factory_.Create(PrecType, matrix_.get(), OverlapLevel));
    assert(ifpackInverse_.get()!=0);
    
    // specify parameters for ILU
    list.set("fact: drop tolerance", 1e-9);  // these should be input parameters from Tramonto
    list.set("fact: level-of-fill", 2);
    //list.set("relaxation: type", "Gauss-Seidel");
    
    // sets the parameters
    IFPACK_CHK_ERR(ifpackInverse_->SetParameters(list));
    
    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    IFPACK_CHK_ERR(ifpackInverse_->Initialize());
  }
  // Builds the preconditioners, by looking for the values of 
  // the matrix.
  IFPACK_CHK_ERR(ifpackInverse_->Compute());

  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_Poisson_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  return(ifpackInverse_->ApplyInverse(X, Y));
}
//==============================================================================
int dft_Poisson_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  return(matrix_->Apply(X, Y));
}
