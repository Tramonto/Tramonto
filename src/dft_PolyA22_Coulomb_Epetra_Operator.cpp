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
  : dft_PolyA22_Epetra_Operator(cmsMap, densityMap, block2Map, parameterList),
    poissonMap_(poissonMap),
    cmsDensMap_(cmsDensMap)
{
  // Construct and label matrices
  poissonOnPoissonMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, poissonMap, 0));
  cmsOnPoissonMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, cmsMap, 0));
  poissonOnDensityMatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, poissonMap, 0));
  Label_ = "dft_PolyA22_Coulomb_Epetra_Operator";
  cmsOnDensityMatrix_->SetLabel("PolyA22Coulomb::cmsOnDensityMatrix");
  cmsOnCmsMatrix_->SetLabel("PolyA22Coulomb::cmsOnCmsMatrix");
  poissonOnPoissonMatrix_->SetLabel("PolyA22Coulomb::poissonOnPoissonMatrix");
  cmsOnPoissonMatrix_->SetLabel("PolyA22Coulomb::cmsOnPoissonMatrix");
  poissonOnDensityMatrix_->SetLabel("PolyA22Coulomb::poissonOnDensityMatrix");

  // Init Ifpack preconditioner for cmsOnCmsMatrix_
  IFPrec = NULL;
  IFPrecType = "ILU"; // incomplete LU
  IFOverlapLevel = 0; // must be >= 0. This value ignored if Comm.NumProc() == 1
  IFList_.set("fact: drop tolerance", 1e-9);
  IFList_.set("fact: level-of-fill", 4);
  // the combine mode is on the following:
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  // Their meaning is as defined in file Epetra_CombineMode.h
  IFList_.set("schwarz: combine mode", "Add");

  // Init ML
  MLPrec = NULL;
  ML_Epetra::SetDefaults("SA",MLList_);
  MLList_.set("ML output", 0);
  // If running TestSmoothers() in ApplyInverse(), uncomment these
  // MLList_.set("test: IFPACK", false);
  // MLList_.set("test: ML self smoother", false);
  MLList_.set("smoother: sweeps",2);
  MLList_.set("smoother: type","MLS");
  MLList_.set("coarse: sweeps", 6);
  MLList_.set("coarse: type", "MLS");
  MLList_.set("coarse: MLS polynomial order", 3); //3 is default
  // ToDo: Turn on repartitioning for better perfomance at extreme scales
  // Need to be able to pass in node positions for partitioning
  //MLList_.set("repartition: enable",1);
  //MLList_.set("repartition: Zoltan dimensions",3);

  // Use this to make ML a direct solve
  // int MaxLevels = 1;
  // MLList_.set("ML output", 10);
  // MLList_.set("max levels", MaxLevels);
  // MLList_.set("coarse: type", "Amesos-KLU");
  
  // Or use Gauss-Seidel 
  /*  
  int MaxLevels = 10;
  int sweeps = 1;
  double omega = 0.67;
  char parameter[80];
  MLList_.set("max levels", MaxLevels);
  for (int ilevel = 0; ilevel < MaxLevels; ilevel++) {
    sprintf(parameter, "smoother: type (level %d)", ilevel);
    MLList_.set(parameter, "Gauss-Seidel");
    sprintf(parameter, "smoother: damping (level %d)", ilevel);
    MLList_.set(parameter, omega);
    sprintf(parameter, "smoother: sweeps (level %d)", ilevel);
    MLList_.set(parameter, sweeps);
  }
  
  MLList_.set("coarse: sweeps", 6);
  MLList_.set("coarse: damping parameter", 0.67); //0.67 is default
  MLList_.set("coarse: type", "Gauss-Seidel");
  */

}
//==============================================================================
dft_PolyA22_Coulomb_Epetra_Operator::~dft_PolyA22_Coulomb_Epetra_Operator() {
  if (MLPrec != 0) delete MLPrec;
  if (IFPrec != 0) delete IFPrec;
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
    poissonOnPoissonMatrix_->PutScalar(0.0);
    cmsOnPoissonMatrix_->PutScalar(0.0);
    poissonOnDensityMatrix_->PutScalar(0.0);
  }
  return(0);
}
//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(int rowGID, int colGID, double value, int blockColFlag) {
  // if poisson then blockColFlag = 0
  // if density then blockColFlag = 1
  // if cms then blockColFlag = 2

  /* The poissonMatrix_, poissonOnDensityMatrix_, cmsOnPoissonMatrix_, and cmsOnDensityMatrix_ values do not change between iterations */  

  if (poissonMap_.MyGID(rowGID)) { // Insert into poissonOnPoissonMatrix or poissonOnDensityMatrix
    if ( blockColFlag == 0 ) { // Insert into poissonOnPoissonMatrix
      if (firstTime_) {
        if (rowGID != curRow_) {
          insertRow();
          curRow_ = rowGID;
        }
        curRowValuesPoissonOnPoisson_[colGID] += value;
      }
      else
        poissonOnPoissonMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
    }
    else if ( blockColFlag == 1 ) { // Insert into poissonOnDensityMatrix
      if (firstTime_) {
        if (rowGID != curRow_) {
          insertRow();
          curRow_ = rowGID;
        }
        curRowValuesPoissonOnDensity_[colGID] += value;
      }
      else
        poissonOnDensityMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(): Invalid argument -- row in poissonMap, but blockColFlag not set for Poisson or density equations.");
      TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } //end if poissonMap_.MyGID(rowGID)
  else if (cmsMap_.MyGID(rowGID)) { // Insert into cmsOnPoissonMatrix or cmsOnCmsMatrix or cmsOnDensityMatrix
    if ( blockColFlag == 0 ) { // Insert into cmsOnPoissonMatrix
      if (firstTime_) {
        if (rowGID != curRow_) {
          insertRow();
          curRow_ = rowGID;
        }
        curRowValuesCmsOnPoisson_[colGID] += value;
      }
      else
        cmsOnPoissonMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
    }
    else if ( blockColFlag == 2 ) { // Insert into cmsOnCmsMatrix
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
        cmsOnDensityMatrix_->SumIntoGlobalValues(rowGID, 1, &value, &colGID);
      }
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(): Invalid argument -- row in cmsMap, but blockColFlag not set for Poisson,density, or cms equations.");
      TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into cmsOnPoisson or cmsOnCmsMatrix or cmsOnDensityMatrix
  else if (densityMap_.MyGID(rowGID)) { // Insert into densityOnDensityMatrix or densityOnCmsMatrix
    if ( blockColFlag == 1 ) { // Insert into densityOnDensityMatrix
      if (rowGID!=colGID) {
        char err_msg[200];
        sprintf(err_msg,"PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(): Invalid argument -- Inserting non-diagonal element into densityOnDensity matrix.");
        TEST_FOR_EXCEPT_MSG(1, err_msg); // Confirm that this is a diagonal value
      }
      (*densityOnDensityMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
    }
    else if ( blockColFlag == 2) { // Insert into densityOnCmsMatrix
      //TEST_FOR_EXCEPT(densityMap_.LID(rowGID)!=cmsMap_.LID(colGID)); // Confirm that this is a diagonal value
      if (densityMap_.LID(rowGID)!=cmsMap_.LID(colGID)) {
        char err_msg[200];
        sprintf(err_msg,"PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(): Invalid argument -- Inserting non-diagonal element into densityOnCms matrix.");
        TEST_FOR_EXCEPT_MSG(1, err_msg); // Confirm that this is a diagonal value
      }
      (*densityOnCmsMatrix_)[densityMap_.LID(rowGID)] += value; // Storing this density block in a vector since it is diagonal
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(): Invalid argument -- row in densityMap, but blockColFlag not set for cms or density equations.");
      TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into densityOnDensityMatrix or densityOnCmsMatrix
  else { // Problem! rowGID not in cmsMap or densityMap or poissonMap
    char err_msg[200];
    sprintf(err_msg,"PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(): rowGID=%i not in cmsMap,densityMap, or poissonMap.",rowGID);
    TEST_FOR_EXCEPT_MSG(1, err_msg);
  }
  
  return(0);
}
//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::insertRow() {

  // Fill matrix rows
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
  if (!curRowValuesPoissonOnPoisson_.empty()) {
    int numEntriesPoissonOnPoisson = curRowValuesPoissonOnPoisson_.size();
    if (numEntriesPoissonOnPoisson>indicesPoissonOnPoisson_.Length()) {
      indicesPoissonOnPoisson_.Resize(numEntriesPoissonOnPoisson);
      valuesPoissonOnPoisson_.Resize(numEntriesPoissonOnPoisson);
    }
    int i=0;
    std::map<int, double>::iterator pos;
    for (pos = curRowValuesPoissonOnPoisson_.begin(); pos != curRowValuesPoissonOnPoisson_.end(); ++pos) {
      indicesPoissonOnPoisson_[i] = pos->first;
      valuesPoissonOnPoisson_[i++] = pos->second;
    }
    poissonOnPoissonMatrix_->InsertGlobalValues(curRow_, numEntriesPoissonOnPoisson, valuesPoissonOnPoisson_.Values(), indicesPoissonOnPoisson_.Values());
  }
  if (!curRowValuesCmsOnPoisson_.empty()) {
    int numEntriesCmsOnPoisson = curRowValuesCmsOnPoisson_.size();
    if (numEntriesCmsOnPoisson>indicesCmsOnPoisson_.Length()) {
      indicesCmsOnPoisson_.Resize(numEntriesCmsOnPoisson);
      valuesCmsOnPoisson_.Resize(numEntriesCmsOnPoisson);
    }
    int i=0;
    std::map<int, double>::iterator pos;
    for (pos = curRowValuesCmsOnPoisson_.begin(); pos != curRowValuesCmsOnPoisson_.end(); ++pos) {
      indicesCmsOnPoisson_[i] = pos->first;
      valuesCmsOnPoisson_[i++] = pos->second;
    }
    cmsOnPoissonMatrix_->InsertGlobalValues(curRow_, numEntriesCmsOnPoisson, valuesCmsOnPoisson_.Values(), indicesCmsOnPoisson_.Values());
  }
  if (!curRowValuesPoissonOnDensity_.empty()) {
    int numEntriesPoissonOnDensity = curRowValuesPoissonOnDensity_.size();
    if (numEntriesPoissonOnDensity>indicesPoissonOnDensity_.Length()) {
      indicesPoissonOnDensity_.Resize(numEntriesPoissonOnDensity);
      valuesPoissonOnDensity_.Resize(numEntriesPoissonOnDensity);
    }
    int i=0;
    std::map<int, double>::iterator pos;
    for (pos = curRowValuesPoissonOnDensity_.begin(); pos != curRowValuesPoissonOnDensity_.end(); ++pos) {
      indicesPoissonOnDensity_[i] = pos->first;
      valuesPoissonOnDensity_[i++] = pos->second;
    }
    poissonOnDensityMatrix_->InsertGlobalValues(curRow_, numEntriesPoissonOnDensity, valuesPoissonOnDensity_.Values(), indicesPoissonOnDensity_.Values());
  }

  curRowValuesCmsOnDensity_.clear();
  curRowValuesCmsOnCms_.clear();
  curRowValuesPoissonOnPoisson_.clear();
  curRowValuesCmsOnPoisson_.clear();
  curRowValuesPoissonOnDensity_.clear();
  return(0);

}

//=============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::finalizeProblemValues() {
  
  if (isLinearProblemSet_) return(0); // nothing to do

  insertRow(); // Dump any remaining entries
  cmsOnCmsMatrix_->FillComplete();
  cmsOnCmsMatrix_->OptimizeStorage();
  
  if (!isFLinear_) {
    cmsOnDensityMatrix_->FillComplete(densityMap_,cmsMap_);
    cmsOnDensityMatrix_->OptimizeStorage();
  }
  
  //only need to do the following the firstTime_
  poissonOnPoissonMatrix_->FillComplete();
  poissonOnPoissonMatrix_->OptimizeStorage();
  cmsOnPoissonMatrix_->FillComplete(poissonMap_, cmsMap_);
  cmsOnPoissonMatrix_->OptimizeStorage();
  poissonOnDensityMatrix_->FillComplete(densityMap_, poissonMap_);
  poissonOnDensityMatrix_->OptimizeStorage();

  if (firstTime_) {
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*poissonOnPoissonMatrix_, MLList_, true);
  }

  if (firstTime_) {
    // allocates an IFPACK factory. No data is associated with this object (only method Create()).
    Ifpack Factory;

    // create the preconditioner. For valid PrecType values, please check the documentation
    IFPrec = Factory.Create(IFPrecType, &(*cmsOnCmsMatrix_), IFOverlapLevel);
    TEST_FOR_EXCEPT(IFPrec == NULL);

    // set the parameters
    IFPACK_CHK_ERR(IFPrec->SetParameters(IFList_));

    // initialize the preconditioner using only filled matrix structure
    // Matrix must have been FillComplete()'d
    IFPACK_CHK_ERR(IFPrec->Initialize());

    // Build the preconditioner using filled matrix values
    IFPACK_CHK_ERR(IFPrec->Compute());
  }

/*
  cout << endl;
  cout << " Number of equations in poissonOnPoisson block = " << poissonOnPoissonMatrix_->NumGlobalRows() << endl;
  cout << " Number of nonzeros in poissonOnPoisson block = " << poissonOnPoissonMatrix_->NumGlobalNonzeros() << endl;
  cout << " Frobenius Norm of poissonOnPoisson block    = " << poissonOnPoissonMatrix_->NormFrobenius() << endl;
  cout << " Average Nonzeros per row of poissonOnPoisson block   = " << ((double)poissonOnPoissonMatrix_->NumGlobalNonzeros())/((double)poissonOnPoissonMatrix_->NumGlobalRows()) << endl;

  cout << endl;
  cout << " Number of equations in poissonOnDensity block = " << poissonOnDensityMatrix_->NumGlobalRows() << endl;
  cout << " Number of nonzeros in poissonOnDensity block = " << poissonOnDensityMatrix_->NumGlobalNonzeros() << endl;
  cout << " Frobenius Norm of poissonOnDensity block    = " << poissonOnDensityMatrix_->NormFrobenius() << endl;
  cout << " Average Nonzeros per row of poissonOnDensity block   = " << ((double)poissonOnDensityMatrix_->NumGlobalNonzeros())/((double)poissonOnDensityMatrix_->NumGlobalRows()) << endl;

  cout << endl;
  cout << " Number of equations in cmsOnPoisson block = " << cmsOnPoissonMatrix_->NumGlobalRows() << endl;
  cout << " Number of nonzeros in cmsOnPoisson block = " << cmsOnPoissonMatrix_->NumGlobalNonzeros() << endl;
  cout << " Frobenius Norm of cmsOnPoisson block    = " << cmsOnPoissonMatrix_->NormFrobenius() << endl;
  cout << " Average Nonzeros per row of cmsOnPoisson block   = " << ((double)cmsOnPoissonMatrix_->NumGlobalNonzeros())/((double)cmsOnPoissonMatrix_->NumGlobalRows()) << endl;

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
int dft_PolyA22_Coulomb_Epetra_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  // If F is in SW (F_location_ == 0):
  // The true A22 block is of the form:

  // |  PP      PD      0    |
  // |  0       DD      DC   |
  // |  CP      CD      CC   |
  
  // where 
  // PP is Poisson on Poisson (general)
  // PD is Poisson on Density (general)
  // DD is Density on Density (diagonal),
  // DC is Density on Cms (diagonal),
  // CP is Cms on Poisson (general)
  // CD is Cms on Density (general, also called the "F matrix")
  // CC is Cms on Cms (general).
  //
  // We will approximate A22 with:
  
  // |  PP      PD       0  |
  // |   0      DD       0  |
  // |  CP      CD      CC  |
  
  // replacing Ddc with a zero matrix for the ApplyInverse method only.

  // Our algorithm is then:
  // Y1 = DD \ X1
  // Y0 = PP \ (X0 - PD*Y1)
  // Y2 = CC \ (X2 - CP*Y0 - CD*Y1)  

  // where inv(DD) is approximated by an ML-generated preconditioner
  // and inv(CC) in approximated using an IFPACK generated preconditioner (currently ILUT)

  // A similar algorithm is found when F is in the NE quadrant

  TEST_FOR_EXCEPT(!X.Map().SameAs(OperatorDomainMap())); 
  TEST_FOR_EXCEPT(!Y.Map().SameAs(OperatorRangeMap()));
  TEST_FOR_EXCEPT(Y.NumVectors()!=X.NumVectors());
  int NumVectors = Y.NumVectors();
  int numCmsElements = cmsMap_.NumMyElements();
  int numDensityElements = densityMap_.NumMyElements(); // == numCmsElements
  int numPoissonElements = poissonMap_.NumMyElements();

  double ** X0ptr;
  double ** Y0ptr;
  double ** X1ptr = new double *[NumVectors];
  double ** Y1ptr = new double *[NumVectors];
  double ** X2ptr = new double *[NumVectors];
  double ** Y2ptr = new double *[NumVectors];

  X.ExtractView(&X0ptr); // Get array of pointers to columns of X
  Y.ExtractView(&Y0ptr); // Get array of pointers to columns of Y
  if (F_location_ == 1) { // F in NE
    for (int i=0; i<NumVectors; i++) {
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      X2ptr[i] = X1ptr[i]+numCmsElements;
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numCmsElements;
    }

    // Hook up pointers into X, Y block vectors
    Epetra_MultiVector X0(View, poissonMap_, X0ptr, NumVectors); // Start X0 to view the first numPoisson elements of X
    Epetra_MultiVector X1(View, cmsMap_, X1ptr, NumVectors); // Start X1 to view middle numCms elements of X
    Epetra_MultiVector X2(View, densityMap_, X2ptr, NumVectors); // Start X2 to view last numDensity elements of X
    Epetra_MultiVector Y0(View, poissonMap_, Y0ptr, NumVectors); // Y0 is a view of the first numPoisson elements of Y
    Epetra_MultiVector Y1(View, cmsMap_, Y1ptr, NumVectors); // Y1 is a view of the middle numCms elements of Y
    Epetra_MultiVector Y2(View, densityMap_, Y2ptr, NumVectors); // Y2 is a view of the last numDensity elements of Y

    // Temporary vectors needed for intermediate results
    Epetra_MultiVector Y0tmp(Y0);
    Epetra_MultiVector Y1tmp1(Y2);
    Epetra_MultiVector Y1tmp2(Y2);

    // Third block row: Y2 = DD\X2
    TEST_FOR_EXCEPT(Y2.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X2, 0.0));
    // First block row: Y0 = PP \ (X0 - PD*Y2);
    TEST_FOR_EXCEPT(poissonOnDensityMatrix_->Apply(Y2, Y0tmp));
    TEST_FOR_EXCEPT(Y0tmp.Update( 1.0, X0, -1.0 ));
    TEST_FOR_EXCEPT(MLPrec->ApplyInverse(Y0tmp, Y0));
    // Third block row: Y1 = CC \ (X1 - CP*Y0 - CD*Y2)
    TEST_FOR_EXCEPT(cmsOnPoissonMatrix_->Apply(Y0, Y1tmp1));
    TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(Y2, Y1tmp2));
    TEST_FOR_EXCEPT( Y1tmp1.Update( 1.0, X1, -1.0, Y1tmp2, -1.0 ) );
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    //Epetra_Vector cmsOnCmsDiag(cmsMap_);
    //cmsOnCmsMatrix_->ExtractDiagonalCopy(cmsOnCmsDiag);
    //TEST_FOR_EXCEPT(Y1.ReciprocalMultiply(1.0, cmsOnCmsDiag, Y1tmp1, 0.0));
    TEST_FOR_EXCEPT(IFPrec->ApplyInverse(Y1tmp1,Y1));

  }
  else { // F in SW
    for (int i=0; i<NumVectors; i++) {
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      X2ptr[i] = X1ptr[i]+numDensityElements;
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numDensityElements;
    }

    // Hook up pointers into X, Y block vectors
    Epetra_MultiVector X0(View, poissonMap_, X0ptr, NumVectors); // Start X0 to view the first numPoisson elements of X
    Epetra_MultiVector X1(View, densityMap_, X1ptr, NumVectors); // Start X1 to view middle numDensity elements of X
    Epetra_MultiVector X2(View, cmsMap_, X2ptr, NumVectors); // Start X2 to view last numCms elements of X
    Epetra_MultiVector Y0(View, poissonMap_, Y0ptr, NumVectors); // Y0 is a view of the first numPoisson elements of Y
    Epetra_MultiVector Y1(View, densityMap_, Y1ptr, NumVectors); // Y1 is a view of the middle numDensity elements of Y
    Epetra_MultiVector Y2(View, cmsMap_, Y2ptr, NumVectors); // Y2 is a view of the last numCms elements of Y

    // Temporary vectors needed for intermediate results
    Epetra_MultiVector Y0tmp(Y0);
    Epetra_MultiVector Y2tmp1(Y2);
    Epetra_MultiVector Y2tmp2(Y2);

    // Second block row: Y1 = DD\X1
    TEST_FOR_EXCEPT(Y1.ReciprocalMultiply(1.0, *densityOnDensityMatrix_, X1, 0.0));
    // First block row: Y0 = PP \ (X0 - PD*Y1);
    TEST_FOR_EXCEPT(poissonOnDensityMatrix_->Apply(Y1, Y0tmp));
    TEST_FOR_EXCEPT(Y0tmp.Update( 1.0, X0, -1.0 ));
    TEST_FOR_EXCEPT(MLPrec->ApplyInverse(Y0tmp, Y0));
    // Third block row: Y2 = CC \ (X2 - CP*Y0 - CD*Y1)
    TEST_FOR_EXCEPT(cmsOnPoissonMatrix_->Apply(Y0, Y2tmp1));
    TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(Y1, Y2tmp2));
    TEST_FOR_EXCEPT( Y2tmp1.Update( 1.0, X2, -1.0, Y2tmp2, -1.0 ) );
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    //Epetra_Vector cmsOnCmsDiag(cmsMap_);
    //cmsOnCmsMatrix_->ExtractDiagonalCopy(cmsOnCmsDiag);
    //TEST_FOR_EXCEPT(Y2.ReciprocalMultiply(1.0, cmsOnCmsDiag, Y2tmp1, 0.0));
    TEST_FOR_EXCEPT(IFPrec->ApplyInverse(Y2tmp1,Y2));

 }

 //To export ML preconditioner
 /* 
   Epetra_MultiVector eye(poissonMap_, poissonMap_.NumGlobalElements(), true);
   Epetra_MultiVector res(poissonMap_, poissonMap_.NumGlobalElements(), true);
   for (int i = 0; i<poissonMap_.NumGlobalElements(); i++) {
   eye.ReplaceGlobalValue(i+poissonMap_.MinAllGID(), i, 1.0);
   }

   MLPrec->ApplyInverse(eye, res);
   EpetraExt::MultiVectorToMatrixMarketFile("prec.dat", res, "", "");
   abort();
 */

  delete [] Y1ptr;
  delete [] X1ptr;
  delete [] Y2ptr;
  delete [] X2ptr;
    
  return(0);
}
//==============================================================================
int dft_PolyA22_Coulomb_Epetra_Operator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {

  // If F is in SW (F_location_ == 0):
  // The A22 block is of the form:

  // |  PP      PD      0    |
  // |  0       DD      DC   |
  // |  CP      CD      CC   |

  // where
  // PP is Poisson on Poisson (general)
  // PD is Poisson on Density (general)
  // DD is Density on Density (diagonal),
  // DC is Density on Cms (diagonal),
  // CP is Cms on Poisson (general)
  // CD is Cms on Density (general, also called the "F matrix")
  // CC is Cms on Cms (general).

  // If F is in NE (F_location_ == 1):
  // The A22 block is of the form:

  // |  PP       0      PD   |
  // |  CP      CC      CD   |
  // |   0      DC      DD   |

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
  if (F_location_ == 1) { // F in NE
    for (int i=0; i<NumVectors; i++) {
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      X2ptr[i] = X1ptr[i]+numCmsElements;
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numCmsElements;
    }
    Epetra_MultiVector X0(View, poissonMap_, X0ptr, NumVectors); // Start X0 to view the first numPoisson elements of X
    Epetra_MultiVector X1(View, cmsMap_, X1ptr, NumVectors); // Start X1 to view middle numCms elements of X
    Epetra_MultiVector X2(View, densityMap_, X2ptr, NumVectors); // Start X2 to view last numDensity elements of X
    Epetra_MultiVector Y0(View, poissonMap_, Y0ptr, NumVectors); // Y0 is a view of the first numPoisson elements of Y
    Epetra_MultiVector Y1(View, cmsMap_, Y1ptr, NumVectors); // Y1 is a view of the middle numCms elements of Y
    Epetra_MultiVector Y2(View, densityMap_, Y2ptr, NumVectors); // Y2 is a view of the last numDensity elements of Y

    // First block row
    Epetra_MultiVector Y0tmp(Y0);
    TEST_FOR_EXCEPT(poissonOnPoissonMatrix_->Apply(X0, Y0));
    TEST_FOR_EXCEPT(poissonOnDensityMatrix_->Apply(X2, Y0tmp));
    TEST_FOR_EXCEPT(Y0.Update(1.0,Y0tmp,1.0));
    // Second block row
    TEST_FOR_EXCEPT(cmsOnPoissonMatrix_->Apply(X0, Y1));
    Epetra_MultiVector Y1tmp1(Y1);
    TEST_FOR_EXCEPT(cmsOnCmsMatrix_->Apply(X1, Y1tmp1));
    Epetra_MultiVector Y1tmp2(Y1);
    TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(X2, Y1tmp2));
    TEST_FOR_EXCEPT(Y1.Update(1.0, Y1tmp1, 1.0));
    TEST_FOR_EXCEPT(Y1.Update(1.0, Y1tmp2, 1.0));
    // Third block row
    TEST_FOR_EXCEPT(Y2.Multiply(1.0, *densityOnCmsMatrix_, X1, 0.0));
    TEST_FOR_EXCEPT(Y2.Multiply(1.0, *densityOnDensityMatrix_, X2, 1.0));
  }
  else { // F in SW
    for (int i=0; i<NumVectors; i++) {
      X1ptr[i] = X0ptr[i]+numPoissonElements;
      X2ptr[i] = X1ptr[i]+numDensityElements;
      Y1ptr[i] = Y0ptr[i]+numPoissonElements;
      Y2ptr[i] = Y1ptr[i]+numDensityElements;
    }

    Epetra_MultiVector X0(View, poissonMap_, X0ptr, NumVectors); // Start X0 to view the first numPoisson elements of X
    Epetra_MultiVector X1(View, densityMap_, X1ptr, NumVectors); // Start X1 to view middle numDensity elements of X
    Epetra_MultiVector X2(View, cmsMap_, X2ptr, NumVectors); // Start X2 to view last numCms elements of X
    Epetra_MultiVector Y0(View, poissonMap_, Y0ptr, NumVectors); // Y0 is a view of the first numPoisson elements of Y
    Epetra_MultiVector Y1(View, densityMap_, Y1ptr, NumVectors); // Y1 is a view of the middle numDensity elements of Y
    Epetra_MultiVector Y2(View, cmsMap_, Y2ptr, NumVectors); // Y2 is a view of the last numCms elements of Y

    // First block row
    Epetra_MultiVector Y0tmp(Y0);
    TEST_FOR_EXCEPT(poissonOnPoissonMatrix_->Apply(X0, Y0));
    TEST_FOR_EXCEPT(poissonOnDensityMatrix_->Apply(X1, Y0tmp));
    TEST_FOR_EXCEPT(Y0.Update(1.0,Y0tmp,1.0));
    // Second block row
    TEST_FOR_EXCEPT(Y1.Multiply(1.0, *densityOnDensityMatrix_, X1, 0.0));
    TEST_FOR_EXCEPT(Y1.Multiply(1.0, *densityOnCmsMatrix_, X2, 1.0));
    // Third block row
    TEST_FOR_EXCEPT(cmsOnPoissonMatrix_->Apply(X0, Y2));
    Epetra_MultiVector Y2tmp1(Y2);
    TEST_FOR_EXCEPT(cmsOnDensityMatrix_->Apply(X1, Y2tmp1));
    Epetra_MultiVector Y2tmp2(Y2);
    TEST_FOR_EXCEPT(cmsOnCmsMatrix_->Apply(X2, Y2tmp2));
    TEST_FOR_EXCEPT(Y2.Update(1.0, Y2tmp1, 1.0));
    TEST_FOR_EXCEPT(Y2.Update(1.0, Y2tmp2, 1.0));
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
  Epetra_Vector bb(b);
  x.Random(); // Fill x with random numbers

  Epetra_Vector x0(View, poissonMap_, x.Values());
  Epetra_Vector b0(View, poissonMap_, b.Values());
  Epetra_Vector x1(View, densityMap_, x.Values()+poissonMap_.NumMyElements()); // Start x1 to view middle numDensity elements of x
  Epetra_Vector b1(View, densityMap_, b.Values()+poissonMap_.NumMyElements()); // Start b1 to view middle numDensity elements of b
  Epetra_Vector x2(View, cmsMap_, x.Values()+poissonMap_.NumMyElements()+densityMap_.NumMyElements()); //Start x2 to view last numCms elements of x
  Epetra_Vector b2(View, cmsMap_, b.Values()+poissonMap_.NumMyElements()+densityMap_.NumMyElements()); // Start b2 to view last numCms elements of b

  // The poisson-on-poisson matrix is singular. Make x0 orthogonal to the kernel.
  Epetra_Vector ones(b0);
  double alpha;
  ones.PutScalar(1.0);
  ones.Norm2(&alpha);
  alpha = 1.0 / alpha;
  ones.Scale( alpha );
  ones.Dot( x0, &alpha );
  alpha = -1.0*alpha;
  x0.Update  ( alpha, ones, 1.0 ); 

/* // Code to test the ML preconditioner
  poissonOnPoissonMatrix_->Apply(x0, b0);
  Epetra_Vector myb(b0);
  MLPrec->ApplyInverse(b0, myb);
  myb.Update(-1.0, x0, 1.0); // Should be zero
  double myresid = 0.0;
  double myxnorm = 0.0;
  myb.Norm2(&myresid);
  x0.Norm2(&myxnorm);
  std::cout << "norm(myb) = " << myresid << endl;
  std::cout << "norm(x0) = " << myxnorm << endl;
*/

  Apply(x, b); // Forward operation
  
  // Inverse if not exact, so we must modify b first:
  if (F_location_ == 1) { // F in NE
    b2.Multiply(-1.0, *densityOnCmsMatrix_, x1, 1.0);
  }
  else { // F in SW
    b1.Multiply(-1.0, *densityOnCmsMatrix_, x2, 1.0);
  }

  ApplyInverse(b, bb); // Reverse operation
  bb.Update(-1.0, x, 1.0); // Should be zero
  double resid = 0.0;
  bb.Norm2(&resid);

  // Warning: poisson-on-poisson matrix is singular; expect low accuracy from
  // ML preconditioner
  if (verbose) 
    std::cout << "A22 self-check residual = " << resid << endl;

  if (resid > 1.0E-12) return(-1); // Bad residual
  return(0);
}
