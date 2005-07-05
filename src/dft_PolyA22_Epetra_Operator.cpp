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

#include "dft_PolyA22_Epetra_Operator.hpp"
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
dft_PolyA22_Epetra_Operator::dft_PolyA22_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & block2Map) 
  : cmsMap_(cmsMap),
    densityMap_(densityMap),
    block2Map_(block2Map),
    cmsOnDensityMatrix_(Epetra_CrsMatrix(Copy, densityMap, 0)),
    densityOnCmsMatrix_(Epetra_Vector(cmsMap)),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    firstTime_(true) {

  Label_ = "dft_PolyA22_Epetra_Operator";
}
//==============================================================================
dft_PolyA22_Epetra_Operator::~dft_PolyA22_Epetra_Operator() {
}
//=============================================================================
int dft_PolyA22_Epetra_Operator::initializeProblemValues() {
  
  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    cmsOnDensityMatrix_.PutScalar(0.0);
    densityOnCmsMatrix_.PutScalar(0.0);
  }
  
  return(0);
}
//=============================================================================
int dft_PolyA22_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value) {

  if (rowGID==colGID) return(0); // diagonals are 1, we don't store them

  if (densityMap_.MyGID(rowGID)) {
    int newRowGID = densityMap_.GID(cmsMap_.LID(rowGID));
    if (firstTime_)
      cmsOnDensityMatrix_.InsertGlobalValues(newRowGID, 1, &value, &colGID);
    else
      cmsOnDensityMatrix_.SumIntoGlobalValues(newRowGID, 1, &value, &colGID);
  }
  else
    densityOnCmsMatrix_[cmsMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal

  return(0);
}
//=============================================================================
int dft_PolyA22_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do

  cmsOnDensityMatrix_.FillComplete();
  cmsOnDensityMatrix_.OptimizeStorage();
  
  /*  std::cout << cmsOnDensityMatrix_<< std::endl;
  */
  if (firstTime_) {
    Teuchos::ParameterList list;
    // create the preconditioner. For valid PrecType values,
    // please check the documentation
    string PrecType = "ILU"; // incomplete LU
    int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
    // it is ignored.
    
    cmsOnDensityInverse_  = Teuchos::rcp(factory_.Create(PrecType, &cmsOnDensityMatrix_, OverlapLevel));
    assert(cmsOnDensityInverse_.get()!=0);
    
    // specify parameters for ILU
    list.set("fact: drop tolerance", 1e-9);  // these should be input parameters from Tramonto
    list.set("fact: level-of-fill", 1);
    
    // sets the parameters
    IFPACK_CHK_ERR(cmsOnDensityInverse_->SetParameters(list));
    
    // initialize the preconditioner. At this point the matrix must
    // have been FillComplete()'d, but actual values are ignored.
    IFPACK_CHK_ERR(cmsOnDensityInverse_->Initialize());
  }
  // Builds the preconditioners, by looking for the values of 
  // the matrix.
  IFPACK_CHK_ERR(cmsOnDensityInverse_->Compute());

  isLinearProblemSet_ = true;
  firstTime_ = false;
  return(0);
}
//==============================================================================
int dft_PolyA22_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numCmsElements = cmsMap_.NumMyElements();

  double ** Y1ptr;
  double ** X1ptr;
  double ** Y2ptr = new double *[NumVectors];
  double ** X2ptr = new double *[NumVectors];

  Y.ExtractView(&Y1ptr); // Get array of pointers to columns of Y
  X.ExtractView(&X1ptr); // Get array of pointers to columns of X
  for (int i=0; i<NumVectors; i++) {
    Y2ptr[i] = Y1ptr[i]+numCmsElements;
    X2ptr[i] = X1ptr[i]+numCmsElements;
  }
  Epetra_MultiVector Y1(View, cmsMap_, Y1ptr, NumVectors); // Y1 is a view of the first numDensity/numCms elements of Y
  Epetra_MultiVector Y2(View, densityMap_, Y2ptr, NumVectors); // Start Y2 to view last numDensity/numCms elements of Y
  Epetra_MultiVector X1(View, densityMap_, X1ptr, NumVectors); // Start X1 to view first numCmsElements elements of X
  Epetra_MultiVector X2(View, cmsMap_, X2ptr, NumVectors); // Start X2 to view last numDensity elements of X

  Y1.ReciprocalMultiply(1.0, densityOnCmsMatrix_, X2, 0.0);
  cmsOnDensityInverse_->ApplyInverse(Y2, X1);
  delete [] Y2ptr;
  delete [] X2ptr;
  return(0);
}
//==============================================================================
int dft_PolyA22_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numCmsElements = cmsMap_.NumMyElements();

  double ** X1ptr;
  double ** Y1ptr;
  double ** X2ptr = new double *[NumVectors];
  double ** Y2ptr = new double *[NumVectors];

  X.ExtractView(&X1ptr); // Get array of pointers to columns of X
  Y.ExtractView(&Y1ptr); // Get array of pointers to columns of Y
  for (int i=0; i<NumVectors; i++) {
    X2ptr[i] = X1ptr[i]+numCmsElements;
    Y2ptr[i] = Y1ptr[i]+numCmsElements;
  }
  Epetra_MultiVector X1a(View, densityMap_, X1ptr, NumVectors); // X1a is a view of the first numDensity/numCms elements of X
  Epetra_MultiVector X2a(View, cmsMap_, X2ptr, NumVectors); // Start X2a to view last numDensity/numCms elements of X
  Epetra_MultiVector X1b(View, cmsMap_, X1ptr, NumVectors); // X1b is a view of the first numDensity/numCms elements of X
  Epetra_MultiVector X2b(View, densityMap_, X2ptr, NumVectors); // Start X2b to view last numDensity/numCms elements of X
  Epetra_MultiVector Y1(View, densityMap_, Y1ptr, NumVectors); // Start Y1 to view first numCmsElements elements of Y
  Epetra_MultiVector Y2(View, cmsMap_, Y2ptr, NumVectors); // Start Y2 to view last numDensity elements of Y

  cmsOnDensityMatrix_.Apply(X2b, Y1);
  Y1.Update(1.0, X1a, 1.0);
  Y2.Multiply(1.0, densityOnCmsMatrix_, X1b, 0.0);
  Y2.Update(1.0, X2a, 1.0);
  delete [] X2ptr;
  delete [] Y2ptr;


  return(0);
}
