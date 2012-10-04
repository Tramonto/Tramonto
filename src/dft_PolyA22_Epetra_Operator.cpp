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
#include "Teuchos_Assert.hpp"

//=============================================================================
dft_PolyA22_Epetra_Operator::dft_PolyA22_Epetra_Operator(const Epetra_Map & cmsMap, const Epetra_Map & densityMap, const Epetra_Map & block2Map, Teuchos::ParameterList * parameterList) 
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

  cmsOnDensityMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cmsMap, 0));
  cmsOnCmsMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cmsMap, 0));
  cmsOnCmsMatrix_->SetLabel("PolyA22::cmsOnCmsMatrix");
  densityOnDensityMatrix_ = Teuchos::rcp(new Epetra_Vector(densityMap));
  densityOnCmsMatrix_ = Teuchos::rcp(new Epetra_Vector(densityMap));
  Label_ = "dft_PolyA22_Epetra_Operator";
  cmsOnDensityMatrix_->SetLabel("PolyA22::cmsOnDensityMatrix");
  F_location_ = Teuchos::getParameter<int>(*parameterList_, "F_location"); //F in NE if F_location = 1, F in SW otherwise

  // Init Ifpack preconditioner for cmsOnCmsMatrix_
  IFPrecType = "ILU"; // incomplete LU
  IFOverlapLevel = 0; // must be >= 0. This value ignored if Comm.NumProc() == 1
  IFList_.set("fact: level-of-fill", 4);

  }
//==============================================================================
dft_PolyA22_Epetra_Operator::~dft_PolyA22_Epetra_Operator() {
}
//=============================================================================
int dft_PolyA22_Epetra_Operator::initializeProblemValues() {

  if (isGraphStructureSet_) return(-1); // Graph structure must be set
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) {
    if (!isFLinear_) cmsOnDensityMatrix_->PutScalar(0.0);
    cmsOnCmsMatrix_->PutScalar(0.0);
    densityOnDensityMatrix_->PutScalar(0.0);
    densityOnCmsMatrix_->PutScalar(0.0);
  }
  
  return(0);
}
//=============================================================================
int dft_PolyA22_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value, int blockColFlag) {
  // if poisson then blockColFlag = 0
  // if density then blockColFlag = 1
  // if cms then blockColFlag = 2

  // Flag to indicate if we've alreay thrown a warning
  static bool offDiagonalDiscarded = false;

  if (cmsMap_.MyGID(rowGID)) { // Insert into cmsOnCmsMatrix or cmsOnDensityMatrix
    if ( blockColFlag == 2 ) { // Insert into cmsOnCmsMatrix
      if (firstTime_) {
	if (rowGID!=curRow_) { 
	  insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
	  curRow_=rowGID;
	}
	curRowValuesCmsOnCms_[colGID] += value;
      }
      else
        cmsOnCmsMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
    }
    else if (blockColFlag == 1) { // Insert into cmsOnDensityMatrix ("F matrix")
      if (firstTime_) {
	if (rowGID!=curRow_) { 
	  insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
	  curRow_=rowGID;
	}
	curRowValuesCmsOnDensity_[colGID] += value;
      }
      else if (!isFLinear_) { 
        //cout<< "row GID = " << rowGID << " value = " << value << " colGID = " << colGID << endl;
	cmsOnDensityMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
      }
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Epetra_Operator::insertMatrixValue(): Invalid argument -- row in cmsMap, but blockColFlag not set for cms or density equations.");
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into cmsOnCmsMatrix or cmsOnDensityMatrix
  else if (densityMap_.MyGID(rowGID)) { // Insert into densityOnDensityMatrix or densityOnCmsMatrix
    if ( blockColFlag == 1 ) { // Insert into densityOnDensityMatrix
      TEUCHOS_TEST_FOR_EXCEPT(rowGID!=colGID); // Confirm that this is a diagonal value
      (*densityOnDensityMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
    }
    else if ( blockColFlag == 2) { // Insert into densityOnCmsMatrix
      //TEUCHOS_TEST_FOR_EXCEPT(densityMap_.LID(rowGID)!=cmsMap_.LID(colGID)); // Confirm that this is a diagonal value
      // The density-on-cms matrix is presumed to be diagonal and stored as a vector. New functionality in the physics breaks this assumption
      // If non-diagonal entries are inserted, discard them and warn the user. This will result in an approximate Jacobian
      if ( densityMap_.LID(rowGID) == cmsMap_.LID(colGID) )
        (*densityOnCmsMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
      else {
       if (!offDiagonalDiscarded) cout << "Warning! Off-diagonal entries detected in density-on-cms block. Jacobian will be inexact.  If code fails try turning off Schur solver." << std::endl;
       offDiagonalDiscarded = true;
      }
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Epetra_Operator::insertMatrixValue(): Invalid argument -- row in densityMap, but blockColFlag not set for cms or density equations.");
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into densityOnDensityMatrix or densityOnCmsMatrix
  else { // Problem! rowGID not in cmsMap or densityMap
    char err_msg[200];
    sprintf(err_msg,"PolyA22_Epetra_Operator::insertMatrixValue(): rowGID=%i not in cmsMap or densityMap.",rowGID);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
  }

  return(0);
}
//=============================================================================
int dft_PolyA22_Epetra_Operator::insertRow() {

  // Fill row of cmsOnCms and cmsOnDensity matrices
  if (!curRowValuesCmsOnDensity_.empty()) {
    int numEntriesCmsOnDensity = curRowValuesCmsOnDensity_.size();
    if (numEntriesCmsOnDensity>indicesCmsOnDensity_.Length()) {
      indicesCmsOnDensity_.Resize(numEntriesCmsOnDensity);
      valuesCmsOnDensity_.Resize(numEntriesCmsOnDensity);
    }
    int i=0;
    std::map<int, double>::iterator pos;
    for (pos = curRowValuesCmsOnDensity_.begin(); pos != curRowValuesCmsOnDensity_.end(); ++pos) {
      indicesCmsOnDensity_[i] = pos->first;
      valuesCmsOnDensity_[i++] = pos->second;
    }
    cmsOnDensityMatrix_->InsertGlobalValues(curRow_, numEntriesCmsOnDensity, valuesCmsOnDensity_.Values(), indicesCmsOnDensity_.Values());
  }
  if (!curRowValuesCmsOnCms_.empty()) {
    int numEntriesCmsOnCms = curRowValuesCmsOnCms_.size();
    if (numEntriesCmsOnCms>indicesCmsOnCms_.Length()) {
      indicesCmsOnCms_.Resize(numEntriesCmsOnCms);
      valuesCmsOnCms_.Resize(numEntriesCmsOnCms);
    }
    int i=0;
    std::map<int, double>::iterator pos;
    for (pos = curRowValuesCmsOnCms_.begin(); pos != curRowValuesCmsOnCms_.end(); ++pos) {
      indicesCmsOnCms_[i] = pos->first;
      valuesCmsOnCms_[i++] = pos->second;
    }
    cmsOnCmsMatrix_->InsertGlobalValues(curRow_, numEntriesCmsOnCms, valuesCmsOnCms_.Values(), indicesCmsOnCms_.Values());
  }
  
  curRowValuesCmsOnDensity_.clear();
  curRowValuesCmsOnCms_.clear();
  return(0);

}
//=============================================================================
int dft_PolyA22_Epetra_Operator::finalizeProblemValues() {
  if (isLinearProblemSet_) return(0); // nothing to do
  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_.NumGlobalElements()==densityMap_.NumGlobalElements(); 

  insertRow(); // Dump any remaining entries
  cmsOnCmsMatrix_->FillComplete();
  cmsOnCmsMatrix_->OptimizeStorage();
  if (!isFLinear_) {
    insertRow(); // Dump any remaining entries
    cmsOnDensityMatrix_->FillComplete(densityMap_, cmsMap_);
    cmsOnDensityMatrix_->OptimizeStorage();
    // cout << " Number of equations in F block = " << cmsOnDensityMatrix_->NumGlobalRows() << endl;
    // cout << " Number of nonzeros in F block = " << cmsOnDensityMatrix_->NumGlobalNonzeros() << endl;
    // cout << " Frobenius Norm of F block    = " << cmsOnDensityMatrix_->NormFrobenius() << endl;
    // cout << " Average Nonzeros per row of F block   = " << cmsOnDensityMatrix_->NumGlobalNonzeros()/cmsOnDensityMatrix_->NumGlobalRows() << endl;
  }

  if (!hasDensityOnCms) { // Confirm that densityOnCmsMatrix is zero
    double normvalue;
    densityOnCmsMatrix_->NormInf(&normvalue);
    TEUCHOS_TEST_FOR_EXCEPT(normvalue!=0.0);
  }

  // IFPack preconditioner for CC block disabled for now
  /*
  if (firstTime_) {
    // allocates an IFPACK factory. No data is associated with this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values, please check the documentation
    IFPrec = Teuchos::rcp( Factory.Create(IFPrecType, &(*cmsOnCmsMatrix_), IFOverlapLevel) );
    TEUCHOS_TEST_FOR_EXCEPT(IFPrec == Teuchos::null);

    // set the parameters
    IFPACK_CHK_ERR(IFPrec->SetParameters(IFList_));

    // initialize the preconditioner using only filled matrix structure
    // Matrix must have been FillComplete()'d
    IFPACK_CHK_ERR(IFPrec->Initialize());

    // Build the preconditioner using filled matrix values
    IFPACK_CHK_ERR(IFPrec->Compute());
  }
  */

/*
  cout << endl;
  cout << " Number of equations in cmsOnDensity block = " << cmsOnDensityMatrix_->NumGlobalRows() << endl;
  cout << " Number of nonzeros in cmsOnDensity block = " << cmsOnDensityMatrix_->NumGlobalNonzeros() << endl;
  cout << " Frobenius Norm of cmsOnDensity block    = " << cmsOnDensityMatrix_->NormFrobenius() << endl;
  cout << " Average Nonzeros per row of cmsOnDensity block   = " << ((double)cmsOnDensityMatrix_->NumGlobalNonzeros())/((double)cmsOnDensityMatrix_->NumGlobalRows()) << endl;

  cout << endl;
  cout << " Number of equations in cmsOnCms block = " << cmsOnCmsMatrix_->NumGlobalRows() << endl;
  cout << " Number of nonzeros in cmsOnCms block = " << cmsOnCmsMatrix_->NumGlobalNonzeros() << endl;
  cout << " Frobenius Norm of cmsOnCms block    = " << cmsOnCmsMatrix_->NormFrobenius() << endl;
  cout << " Average Nonzeros per row of cmsOnCms block   = " << ((double)cmsOnCmsMatrix_->NumGlobalNonzeros())/((double)cmsOnCmsMatrix_->NumGlobalRows()) << endl;

  cout << endl;
  double norm;
  densityOnDensityMatrix_->Norm2( &norm  );
  cout << " Number of equations in densityOnDensity block = " << densityOnDensityMatrix_->GlobalLength() << endl;
  cout << " Number of nonzeros in densityOnDensity block = " << densityOnDensityMatrix_->GlobalLength() << endl;
  cout << " Frobenius Norm of densityOnDensity block    = " << norm << endl;
  cout << " Average Nonzeros per row of DensityOnDensity block   = " << 1.0 << endl;

  cout << endl;
  densityOnCmsMatrix_->Norm2( &norm  );
  cout << " Number of equations in densityOnCms block = " << densityOnCmsMatrix_->GlobalLength() << endl;
  cout << " Number of nonzeros in densityOnCms block = " << densityOnCmsMatrix_->GlobalLength() << endl;
  cout << " Frobenius Norm of densityOnCms block    = " << norm << endl;
  cout << " Average Nonzeros per row of DensityOnCms block   = " << 1.0 << endl;
  cout << endl;
*/

  isLinearProblemSet_ = true;
  firstTime_ = false;

  return(0);
}
//==============================================================================
int dft_PolyA22_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  // If F is in SW (F_location_ == 0):
  // The true A22 block is of the form:

  // | DD      DC   |
  // | CD      CC   |

  // where
  // DD is Density on Density (diagonal),
  // DC is Density on Cms (diagonal),
  // CD is Cms on Density (general, also called the "F matrix")
  // CC is Cms on Cms (general).
  //
  // We will approximate A22 with:

  // | DD       0  |
  // | CD      CC  |

  // replacing DC with a zero matrix for the ApplyInverse method only.

  // Our algorithm is then:
  // Y1 = DD \ X1
  // Y2 = CC \ (X2 - CD*Y1)

  // where inv(DD) is approximated by an ML-generated preconditioner
  // and inv(CC) is approximated using a diagonal preconditioner
  // (Code to use IFPACK generated preconditioner (currently ILU) commented out.)

  // A similar algorithm is found when F is in the NE quadrant

  //double normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in PolyA22 ApplyInverse = " << normvalue << endl;

  TEUCHOS_TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEUCHOS_TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
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
    
    // Temporary vectors needed for intermediate results
    Epetra_MultiVector Y1tmp(Y1);

    // Second block row: Y2 = DD\X2
    TEUCHOS_TEST_FOR_EXCEPT(Y2.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X2, 0.0));
    // First block row: Y1 = CC \ (X1 - CD*Y2)
    TEUCHOS_TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(Y2, Y1tmp));
    TEUCHOS_TEST_FOR_EXCEPT(Y1tmp.Update(1.0, X1, -1.0));
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    Epetra_Vector cmsOnCmsDiag(cmsMap_);
    cmsOnCmsMatrix_->ExtractDiagonalCopy(cmsOnCmsDiag);
    TEUCHOS_TEST_FOR_EXCEPT(Y1.ReciprocalMultiply(1.0, cmsOnCmsDiag, Y1tmp, 0.0));
    // TEUCHOS_TEST_FOR_EXCEPT(IFPrec->ApplyInverse(Y1tmp,Y1));
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

    // Temporary vectors needed for intermediate results
    Epetra_MultiVector Y2tmp(Y2);

    // First block row: Y1 = DD\X1
    TEUCHOS_TEST_FOR_EXCEPT(Y1.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X1, 0.0));
    // Second block row: Y2 = CC \ (X2 - CD*Y1)
    TEUCHOS_TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(Y1, Y2tmp));
    TEUCHOS_TEST_FOR_EXCEPT(Y2tmp.Update(1.0, X2, -1.0));
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    Epetra_Vector cmsOnCmsDiag(cmsMap_);
    cmsOnCmsMatrix_->ExtractDiagonalCopy(cmsOnCmsDiag);
    TEUCHOS_TEST_FOR_EXCEPT(Y2.ReciprocalMultiply(1.0, cmsOnCmsDiag, Y2tmp, 0.0));
    // TEUCHOS_TEST_FOR_EXCEPT(IFPrec->ApplyInverse(Y2tmp,Y2));
  }
  
  delete [] Y2ptr;
  delete [] X2ptr;

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA22 ApplyInverse = " << normvalue << endl;
  return(0);
}
//==============================================================================
int dft_PolyA22_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  // If F is in SW (F_location_ == 0):
  // The A22 block is of the form:

  // |  DD      DC   |
  // |  CD      CC   |

  // where
  // DD is Density on Density (diagonal),
  // DC is Density on Cms (diagonal),
  // CD is Cms on Density (general, also called the "F matrix")
  // CC is Cms on Cms (general).

  // If F is in NE (F_location_ == 1):
  // The A22 block is of the form:

  // |  CC      CD   |
  // |  DC      DD   |

  //double normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in PolyA22 Apply = " << normvalue << endl;

  TEUCHOS_TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
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

    TEUCHOS_TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(X2, Y1));
    Epetra_MultiVector Y1tmp(Y1);
    TEUCHOS_TEST_FOR_EXCEPT(cmsOnCmsMatrix_->Apply(X1, Y1tmp));
    TEUCHOS_TEST_FOR_EXCEPT(Y1.Update(1.0, Y1tmp, 1.0));
    TEUCHOS_TEST_FOR_EXCEPT(Y2.Multiply(1.0, *densityOnCmsMatrix_, X1, 0.0));
    TEUCHOS_TEST_FOR_EXCEPT(Y2.Multiply(1.0, *densityOnDensityMatrix_, X2, 1.0));
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
      TEUCHOS_TEST_FOR_EXCEPT(Y1.Multiply(1.0, *densityOnCmsMatrix_, X2tmp, 0.0)); // Momentarily make X2 compatible with densityMap_
    }
    TEUCHOS_TEST_FOR_EXCEPT(Y1.Multiply(1.0, *densityOnDensityMatrix_, X1, 1.0));
    TEUCHOS_TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(X1, Y2));
    Epetra_MultiVector Y2tmp(Y2);
    TEUCHOS_TEST_FOR_EXCEPT(cmsOnCmsMatrix_->Apply(X2, Y2tmp));
    TEUCHOS_TEST_FOR_EXCEPT(Y2.Update(1.0, Y2tmp, 1.0));
  }
  
  delete [] X2ptr;
  delete [] Y2ptr;

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA22 Apply = " << normvalue << endl;
  return(0);
}
//==============================================================================
int dft_PolyA22_Epetra_Operator::Check(bool verbose) const {

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
