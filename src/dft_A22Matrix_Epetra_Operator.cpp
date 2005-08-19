//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
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

#include "dft_A22Matrix_Epetra_Operator.hpp"
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
dft_A22Matrix_Epetra_Operator::dft_A22Matrix_Epetra_Operator(const Epetra_Map & block2Map) 
  : block2Map_(block2Map),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true) {

  A22Matrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, block2Map, 0));
  Label_ = "dft_A22Matrix_Epetra_Operator";
  A22Matrix_->SetLabel("dft_A22Matrix_Epetra_Operator::A22Matrix");
}
//==============================================================================
dft_A22Matrix_Epetra_Operator::~dft_A22Matrix_Epetra_Operator() {
}
//=============================================================================
int dft_A22Matrix_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_)
    A22Matrix_->PutScalar(0.0);
  
  return(0);
}
//=============================================================================
int dft_A22Matrix_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {
  
  
  if (firstTime_)
    A22Matrix_->InsertGlobalValues(rowGID, 1, &value, &colGID);
  else
    A22Matrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);

  return(0);
}
//=============================================================================
int dft_A22Matrix_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  A22Matrix_->FillComplete();
  A22Matrix_->OptimizeStorage();
  //EpetraExt::RowMatrixToMatrixMarketFile( "A22.dat", A22Matrix_, "CMS and Density blocks", 
  //				  "The 2,2 block of  problem",
  //				  true);
 // abort();
  
  /*  std::cout << A22Matrix_<< std::endl;
  */

  if (firstTime_) {
    Teuchos::ParameterList list;
    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    string PrecType = "ILU"; // incomplete LU
    //string PrecType = "point relaxation"; // Gauss-Seidel
    int OverlapLevel = 0; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.
    
    A22Inverse_  = Teuchos::rcp(factory_.Create(PrecType, A22Matrix_.get(), OverlapLevel));
    assert(A22Inverse_.get()!=0);
    
    // specify parameters for ILU
    list.set("fact: drop tolerance", 1e-9);  // these should be input parameters from Tramonto
    list.set("fact: level-of-fill", 2);
    //list.set("relaxation: type", "Gauss-Seidel");
    
    // sets the parameters
    IFPACK_CHK_ERR(A22Inverse_->SetParameters(list));
    
    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    IFPACK_CHK_ERR(A22Inverse_->Initialize());
  }
  // Builds the preconditioners, by looking for the values of 
  // the matrix.
  IFPACK_CHK_ERR(A22Inverse_->Compute());

  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_A22Matrix_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  return(A22Inverse_->ApplyInverse(X, Y));
}
//==============================================================================
int dft_A22Matrix_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());

  return(A22Matrix_->Apply(X, Y));
}
