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

#include "dft_PolyA22_Coulomb_Tpetra_Operator.hpp"

//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA22_Coulomb_Tpetra_Operator
(const RCP<const MAP> & cmsMap, const RCP<const MAP> & densityMap, 
 const RCP<const MAP> & poissonMap, const RCP<const MAP> & cmsDensMap, 
 const RCP<const MAP> & block2Map,  RCP<ParameterList> parameterList) 
  : dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>
              (cmsMap, densityMap, cmsDensMap, parameterList),
    poissonMap_(poissonMap),
    cmsDensMap_(cmsDensMap),
    block2Map_(block2Map),
    curPoissonRow_(-1),
    curCPRow_(-1),
    curPDRow_(-1)
{
  poissonOnPoissonMatrix_ = rcp(new MAT(poissonMap, 0));    
  cmsOnPoissonMatrix_ = rcp(new MAT(cmsMap, 0));
  poissonOnDensityMatrix_ = rcp(new MAT(poissonMap, 0));
  Label_ = "dft_PolyA22_Coulomb_Tpetra_Operator";
  cmsOnDensityMatrix_->setObjectLabel("PolyA22Coulomb::cmsOnDensityMatrix");
  cmsOnCmsMatrix2_->setObjectLabel("PolyA22Coulomb::cmsOnCmsMatrix");
  poissonOnPoissonMatrix_->setObjectLabel("PolyA22Coulomb::poissonOnPoissonMatrix");
  cmsOnPoissonMatrix_->setObjectLabel("PolyA22Coulomb::cmsOnPoissonMatrix");
  poissonOnDensityMatrix_->setObjectLabel("PolyA22Coulomb::poissonOnDensityMatrix");
  /*
  ML_Epetra::SetDefaults("SA",*MLList_);
  //MLList_->set("ML output", 0);
   If running TestSmoothers() in applyInverse(), uncomment these
     MLList_->set("test: IFPACK", false);
     MLList_->set("test: ML self smoother", false);
  
  int MaxLevels = 10;
  int sweeps = 2; 
  Scalar alpha = 20.0; //30.0 is default
  char parameter[80];
  MLList_->set("max levels", MaxLevels);
  for (LocalOrdinal ilevel = 0; ilevel < MaxLevels; ilevel++) 
  {
    sprintf(parameter, "smoother: type (level %d)", ilevel);
    MLList_->set(parameter, "MLS");
    sprintf(parameter, "smoother: MLS polynomial order (level %d)", ilevel);
    MLList_->set(parameter, 3); //3 is default
    sprintf(parameter, "smoother: MLS alpha (level %d)", ilevel);
    MLList_->set(parameter, alpha);
  } //end for
  
  MLList_->set("coarse: sweeps", 6); 
  MLList_->set("coarse: type", "MLS");
  MLList_->set("coarse: MLS polynomial order", 3); //3 is default
  
  // Or use Gauss-Seidel 
    
  LocalOrdinal MaxLevels = 10;
  LocalOrdinal sweeps = 1;
  Scalar omega = 0.67;
  char parameter[80];
  MLList_->set("max levels", MaxLevels);
  for (LocalOrdinal ilevel = 0; ilevel < MaxLevels; ilevel++) 
  {
    sprintf(parameter, "smoother: type (level %d)", ilevel);
    MLList_->set(parameter, "Gauss-Seidel");
    sprintf(parameter, "smoother: damping (level %d)", ilevel);
    MLList_->set(parameter, omega);
    sprintf(parameter, "smoother: sweeps (level %d)", ilevel);
    MLList_->set(parameter, sweeps);
  } //end for
  
  MLList->set("coarse: sweeps", 6);
  MLList_->set("coarse: damping parameter", 0.67); //0.67 is default
  MLList_->set("coarse: type", "Gauss-Seidel");
  */

} //end constructor
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_PolyA22_Coulomb_Tpetra_Operator
() 
{
  return;
} //end destructor
//============================================================================= 
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
() 
{

  TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n"); 
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) 
  {
    if (!isFLinear_) 
    {
      cmsOnDensityMatrix_->resumeFill();
      cmsOnDensityMatrix_->setAllToScalar(0.0);
    } //end if
    cmsOnCmsMatrix2_->resumeFill();
    cmsOnCmsMatrix2_->setAllToScalar(0.0);
    densityOnDensityMatrix_->putScalar(0.0);
    densityOnCmsMatrix_->putScalar(0.0);
    poissonOnPoissonMatrix_->resumeFill();
    poissonOnPoissonMatrix_->setAllToScalar(0.0);
    cmsOnPoissonMatrix_->resumeFill();
    cmsOnPoissonMatrix_->setAllToScalar(0.0);
    poissonOnDensityMatrix_->resumeFill();
    poissonOnDensityMatrix_->setAllToScalar(0.0);
  } //end if
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value, GlobalOrdinal blockColFlag) 
{
  // if poisson then blockColFlag = 0
  // if density then blockColFlag = 1
  // if cms then blockColFlag = 2

  Array<GlobalOrdinal> cols(1);
  cols[0] = colGID;
  Array<Scalar> vals(1);
  vals[0] = value;

  /* The poissonMatrix_, poissonOnDensityMatrix_, cmsOnPoissonMatrix_, and cmsOnDensityMatrix_ values do not change between iterations */  

  if (poissonMap_->isNodeGlobalElement(rowGID)) { // Insert into poissonOnPoissonMatrix or poissonOnDensityMatrix
    if ( blockColFlag == 0 ) { // Insert into poissonOnPoissonMatrix
      if (firstTime_) {
	if (rowGID != curRow_) {
	  insertRow();
	  curRow_ = rowGID;
	}
	curRowValuesPoissonOnPoisson_[colGID] += value;
      }
      else
	poissonOnPoissonMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
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
      	poissonOnDensityMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- row in poissonMap, but blockColFlag not set for Poisson or density equations.");
      TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } //end if poissonMap_.MyGID(rowGID)
  else if (cmsMap_->isNodeGlobalElement(rowGID)) { // Insert into cmsOnPoissonMatrix or cmsOnCmsMatrix or cmsOnDensityMatrix
    if ( blockColFlag == 0 ) { // Insert into cmsOnPoissonMatrix
      if (firstTime_) {
	if (rowGID != curRow_) {
	  insertRow();
	  curRow_ = rowGID;
	}
	curRowValuesCmsOnPoisson_[colGID] += value;
      }
      else
	cmsOnPoissonMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
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
	cmsOnCmsMatrix2_->sumIntoGlobalValues(rowGID, cols, vals);
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
      	cmsOnDensityMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
      }
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- row in cmsMap, but blockColFlag not set for Poisson,density, or cms equations.");
      TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into cmsOnPoisson or cmsOnCmsMatrix or cmsOnDensityMatrix
  else if (densityMap_->isNodeGlobalElement(rowGID)) { // Insert into densityOnDensityMatrix or densityOnCmsMatrix
    if ( blockColFlag == 1 ) { // Insert into densityOnDensityMatrix
      if (rowGID!=colGID) {
	char err_msg[200];
	sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- Inserting non-diagonal element into densityOnDensity matrix.");
	TEST_FOR_EXCEPT_MSG(1, err_msg); // Confirm that this is a diagonal value
      }
      densityOnDensityMatrix_->sumIntoLocalValue(densityMap_->getLocalElement(rowGID), value);
    }
    else if ( blockColFlag == 2) { // Insert into densityOnCmsMatrix
      if (densityMap_->getLocalElement(rowGID)!=cmsMap_->getLocalElement(colGID)) {
	char err_msg[200];
	sprintf(err_msg,"PolyA22_Coulomb_Epetra_Operator::insertMatrixValue(): Invalid argument -- Inserting non-diagonal element into densityOnCms matrix.");
	TEST_FOR_EXCEPT_MSG(1, err_msg); // Confirm that this is a diagonal value
      }
      densityOnCmsMatrix_->sumIntoLocalValue(densityMap_->getLocalElement(rowGID), value);
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): Invalid argument -- row in densityMap, but blockColFlag not set for cms or density equations.");
      TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into densityOnDensityMatrix or densityOnCmsMatrix
  else { // Problem! rowGID not in cmsMap or densityMap or poissonMap
    char err_msg[200];
    sprintf(err_msg,"PolyA22_Coulomb_Tpetra_Operator::insertMatrixValue(): rowGID=%i not in cmsMap,densityMap, or poissonMap.",rowGID);
    TEST_FOR_EXCEPT_MSG(1, err_msg);
  }

} //end insertMatrixValue
//=============================================================================
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void 
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
() 
{
  // Fill matrix rows
  if (!curRowValuesCmsOnDensity_.empty()) {
    int numEntriesCmsOnDensity = curRowValuesCmsOnDensity_.size();
    if (numEntriesCmsOnDensity>indicesCmsOnDensity_.size()) {
      indicesCmsOnDensity_.resize(numEntriesCmsOnDensity);
      valuesCmsOnDensity_.resize(numEntriesCmsOnDensity);
    }
    LocalOrdinal i=0;
    typename std::map<GlobalOrdinal, Scalar>::iterator pos;
    for (pos = curRowValuesCmsOnDensity_.begin(); pos != curRowValuesCmsOnDensity_.end(); ++pos) {
      indicesCmsOnDensity_[i] = pos->first;
      valuesCmsOnDensity_[i++] = pos->second;
    }
    cmsOnDensityMatrix_->insertGlobalValues(curRow_, indicesCmsOnDensity_, valuesCmsOnDensity_);
  }
  if (!curRowValuesCmsOnCms_.empty()) {
    int numEntriesCmsOnCms = curRowValuesCmsOnCms_.size();
    if (numEntriesCmsOnCms>indicesCmsOnCms_.size()) {
      indicesCmsOnCms_.resize(numEntriesCmsOnCms);
      valuesCmsOnCms_.resize(numEntriesCmsOnCms);
    }
    LocalOrdinal i=0;
    typename std::map<GlobalOrdinal, Scalar>::iterator pos;
    for (pos = curRowValuesCmsOnCms_.begin(); pos != curRowValuesCmsOnCms_.end(); ++pos) {
      indicesCmsOnCms_[i] = pos->first;
      valuesCmsOnCms_[i++] = pos->second;
    }
    cmsOnCmsMatrix2_->insertGlobalValues(curRow_, indicesCmsOnCms_, valuesCmsOnCms_);
  }
  if (!curRowValuesPoissonOnPoisson_.empty()) {
    int numEntriesPoissonOnPoisson = curRowValuesPoissonOnPoisson_.size();
    if (numEntriesPoissonOnPoisson>indicesPoissonOnPoisson_.size()) {
      indicesPoissonOnPoisson_.resize(numEntriesPoissonOnPoisson);
      valuesPoissonOnPoisson_.resize(numEntriesPoissonOnPoisson);
    }
    LocalOrdinal i=0;
    typename std::map<GlobalOrdinal, Scalar>::iterator pos;
    for (pos = curRowValuesPoissonOnPoisson_.begin(); pos != curRowValuesPoissonOnPoisson_.end(); ++pos) {
      indicesPoissonOnPoisson_[i] = pos->first;
      valuesPoissonOnPoisson_[i++] = pos->second;
    }
    poissonOnPoissonMatrix_->insertGlobalValues(curRow_, indicesPoissonOnPoisson_, valuesPoissonOnPoisson_);
  }
  if (!curRowValuesCmsOnPoisson_.empty()) {
    int numEntriesCmsOnPoisson = curRowValuesCmsOnPoisson_.size();
    if (numEntriesCmsOnPoisson>indicesCmsOnPoisson_.size()) {
      indicesCmsOnPoisson_.resize(numEntriesCmsOnPoisson);
      valuesCmsOnPoisson_.resize(numEntriesCmsOnPoisson);
    }
    LocalOrdinal i=0;
    typename std::map<GlobalOrdinal, Scalar>::iterator pos;
    for (pos = curRowValuesCmsOnPoisson_.begin(); pos != curRowValuesCmsOnPoisson_.end(); ++pos) {
      indicesCmsOnPoisson_[i] = pos->first;
      valuesCmsOnPoisson_[i++] = pos->second;
    }
    cmsOnPoissonMatrix_->insertGlobalValues(curRow_, indicesCmsOnPoisson_, valuesCmsOnPoisson_);
  }
  if (!curRowValuesPoissonOnDensity_.empty()) {
    int numEntriesPoissonOnDensity = curRowValuesPoissonOnDensity_.size();
    if (numEntriesPoissonOnDensity>indicesPoissonOnDensity_.size()) {
      indicesPoissonOnDensity_.resize(numEntriesPoissonOnDensity);
      valuesPoissonOnDensity_.resize(numEntriesPoissonOnDensity);
    }
    LocalOrdinal i=0;
    typename std::map<GlobalOrdinal, Scalar>::iterator pos;
    for (pos = curRowValuesPoissonOnDensity_.begin(); pos != curRowValuesPoissonOnDensity_.end(); ++pos) {
      indicesPoissonOnDensity_[i] = pos->first;
      valuesPoissonOnDensity_[i++] = pos->second;
    }
    poissonOnDensityMatrix_->insertGlobalValues(curRow_, indicesPoissonOnDensity_, valuesPoissonOnDensity_);
  }

  indicesCmsOnDensity_.clear();
  valuesCmsOnDensity_.clear();
  indicesCmsOnCms_.clear();
  valuesCmsOnCms_.clear();
  indicesPoissonOnPoisson_.clear();
  valuesPoissonOnPoisson_.clear();
  indicesCmsOnPoisson_.clear();
  valuesCmsOnPoisson_.clear();
  indicesPoissonOnDensity_.clear();
  valuesPoissonOnDensity_.clear();
  curRowValuesCmsOnDensity_.clear();
  curRowValuesCmsOnCms_.clear();
  curRowValuesPoissonOnPoisson_.clear();
  curRowValuesCmsOnPoisson_.clear();
  curRowValuesPoissonOnDensity_.clear();

}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
() 
{
  if (isLinearProblemSet_) 
  {
    return; // nothing to do
  } //end if
  insertRow(); // Dump any remaining entries
  cmsOnCmsMatrix2_->fillComplete();
  
  if (!isFLinear_) {
    cmsOnDensityMatrix_->fillComplete(densityMap_,cmsMap_);
  } //end if<
  
  //only need to do the following the firstTime_
  poissonOnPoissonMatrix_->fillComplete();
  cmsOnPoissonMatrix_->fillComplete(poissonMap_, cmsMap_);
  poissonOnDensityMatrix_->fillComplete(densityMap_, poissonMap_);

  isLinearProblemSet_ = true;
  firstTime_ = false;
} //end finalizeProblemValues
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const 
{
  // If F is in SW (F_location_ == 0):
  // The true A22 block is of the form:

  // |  P       Pd      0    |
  // |  0       Ddd     Ddc  |
  // |  Cp      F       Dcc  |
  
  // where 
  // P is Poisson matrix,
  // Pd is Poisson on Density matrix,
  // Ddd is Density on Density (diagonal),
  // Ddc is Density on Cms (diagonal with small coefficient values),
  // Cp is Cms on Poisson matrix,
  // F is Cms on Density (fairly dense),
  // Dcc is Cms on Cms (diagonal).
  //
  // We will approximate A22 with:
  
  // |  P       Pd      0    |
  // |  0       Ddd     0    |
  // |  Cp      F       Dcc  |
  
  // replacing Ddc with a zero matrix for the applyInverse method only.

  // Our algorithm is then:
  // Y1 = Ddd \ X1
  // Y0 = P \ (X0 - Pd*Y1)
  // Y2 = Dcc \ (X2 - Cp*Y0 - F*Y1)  

  // A similar algorithm is found when F is in the NE quadrant

  TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap())); 
  TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  printf("\n\n\n\ndft_PolyA22_Coulomb_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  size_t NumVectors = Y.getNumVectors();
  size_t numCmsElements = cmsMap_->getNodeNumElements();
  size_t numDensityElements = densityMap_->getNodeNumElements(); // == numCmsElements
  size_t numPoissonElements = poissonMap_->getNodeNumElements();

  RCP<const MV> X0 = X.offsetView(poissonMap_, 0);
  // X0 is a view of the first numPoisson elements of X
  RCP<const MV> X1 = X.offsetView(densityMap_, numPoissonElements);
  // X1 is a view of the middle numDensity/numCms elements of X
  RCP<const MV> X2;
  RCP<MV> Y0 = Y.offsetViewNonConst(poissonMap_, 0);
  // Y0 is a view of the first numPoisson elements of Y
  RCP<MV> Y1 = Y.offsetViewNonConst(densityMap_, numPoissonElements);
  // Y1 is a view of the middle numDensity/numCms elements of Y
  RCP<MV> Y2;

  if (F_location_ == 1)  //F in NE part
  {
    X2 = X.offsetView(densityMap_, numCmsElements);
    Y2 = Y.offsetViewNonConst(densityMap_, numCmsElements);
  } //end if
  else  //F in SW part
  {
    X2 = X.offsetView(densityMap_, numDensityElements);
    Y2 = Y.offsetViewNonConst(densityMap_, numDensityElements);
  } //end else
  // X2 is a view of the last numDensity/numCms elements of X
  // Y2 is a view of the last numDensity/numCms elements of Y

  RCP<MV > Y0tmp = rcp(new MV(*Y0));
  RCP<MV > Y1tmp1 = rcp(new MV(*Y1));
  RCP<MV > Y1tmp2 = rcp(new MV(*Y1));
  RCP<MV > Y2tmp1 = rcp(new MV(*Y2));
  RCP<MV > Y2tmp2 = rcp(new MV(*Y2));

  RCP<VEC> tmp = rcp(new VEC(densityMap_));
  RCP<VEC> tmp2 = rcp(new VEC(densityMap_));
  tmp->reciprocal(*densityOnDensityMatrix_);
  if (F_location_ == 1) 
  {
    // Third block row: Y2 = DD\X2
    Y2->elementWiseMultiply(1.0, *tmp, *X2, 0.0);
    // First block row: Y0 = PP \ (X0 - PD*Y2);
    poissonOnDensityMatrix_->apply(*Y2, *Y0tmp);
    Y0tmp->update(1.0, *X0, -1.0);
    //// equivalent of ML applying inv(PP) ...? ////
    // Third block row: Y1 = CC \ (X1 - CP*Y0 - CD*Y2)
    cmsOnPoissonMatrix_->apply(*Y0, *Y1tmp1);
    cmsOnDensityMatrix_->apply(*Y2, *Y1tmp2);
    Y1tmp1->update(1.0, *X1, -1.0, *Y1tmp2, -1.0);    
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    VEC cmsOnCmsDiag(cmsMap_);
    cmsOnCmsMatrix2_->getLocalDiagCopy(cmsOnCmsDiag);
    tmp2->reciprocal(cmsOnCmsDiag);
    Y1->elementWiseMultiply(1.0, *tmp2, *Y1tmp1, 0.0);

  } //end if
  else 
  {
    // Second block row: Y1 = DD\X1
    Y1->elementWiseMultiply(1.0, *tmp, *X1, 0.0);
    // First block row: Y0 = PP \ (X0 - PD*Y1);
    poissonOnDensityMatrix_->apply(*Y1, *Y0tmp);
    Y0tmp->update( 1.0, *X0, -1.0 );
    //// equivalent of ML applying inv(PP) ...? ////
    // Third block row: Y2 = CC \ (X2 - CP*Y0 - CD*Y1)
    cmsOnPoissonMatrix_->apply(*Y0, *Y2tmp1);
    cmsOnDensityMatrix_->apply(*Y1, *Y2tmp2);
    Y2tmp1->update(1.0, *X2, -1.0, *Y2tmp2, -1.0);
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    VEC cmsOnCmsDiag(cmsMap_);
    cmsOnCmsMatrix2_->getLocalDiagCopy(cmsOnCmsDiag);
    tmp2->reciprocal(cmsOnCmsDiag);
    Y2->elementWiseMultiply(1.0, *tmp2, *Y2tmp1, 0.0);

  } //end else
  /*    
#ifdef SUPPORTS_STRATIMIKOS
  RCP<ThyraMV> thyraY = createMultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(Y0);
  RCP<ThyraOP> thyraOp = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(poissonMatrix_);
  
  RCP<DefaultLinearSolverBuilder> solver = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver->readParameters(out.get());
  
  RCP<ThyraLOWSFactory> lowsFactory = solver->createLinearSolveStrategy("");
  RCP<ThyraLOWS> lows = linearOpWithSolve<Scalar>(*lowsFactory, thyraOp);
  SolveStatus<Scalar> status = lows->solve(Thyra::NOTRANS, *thyraY, thyraY.ptr());
#else
  RCP<LinPROB> problem = rcp(new LinPROB(poissonMatrix_, Y0, Y0));
  TEST_FOR_EXCEPT(problem->setProblem() == false);
  RCP<SolMGR> solver = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem, parameterList_));
  ReturnType ret = solver->solve();
#endif
  */  

} //end applyInverse
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const 
{
  TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
  size_t NumVectors = Y.getNumVectors();
  size_t numCmsElements = cmsMap_->getNodeNumElements();
  size_t numDensityElements = densityMap_->getNodeNumElements();
  size_t numPoissonElements = poissonMap_->getNodeNumElements();

  RCP<const MV> X0 = X.offsetView(poissonMap_, 0);
  // Start X0 to view the first numPoisson elements of X
  RCP<const MV> X1 = X.offsetView(densityMap_, numPoissonElements);
  // Start X1 to view middle numDensityElements/numCms elements of X
  RCP<const MV> X2;
  // Start X2 to view last numDensity/numCms elements of X - was cms
  RCP<MV > Y0 = Y.offsetViewNonConst(poissonMap_, 0);
  // Y0 is a view of the first numPoisson elements of Y
  RCP<MV > Y1 = Y.offsetViewNonConst(densityMap_, numPoissonElements);
  // Y1 is a view of the middle numDensity/numCms elements of Y
  RCP<MV > Y2;
  // Y2 is a view of the last numDensity/numCms elements of Y

  if (F_location_ == 1) 
  {
    X2 = X.offsetView(densityMap_, numCmsElements);
    Y2 = Y.offsetViewNonConst(densityMap_, numCmsElements);
  } //end if
  else 
  {
    X2 = X.offsetView(densityMap_, numDensityElements);
    Y2 = Y.offsetViewNonConst(densityMap_, numDensityElements);
  } //end else

  RCP<MV > Y0tmp = rcp(new MV(*Y0));
  RCP<MV > Y0tmp2 = rcp(new MV(*Y0));
  RCP<MV > Y1tmp1 = rcp(new MV(*Y1));
  RCP<MV > Y1tmp2 = rcp(new MV(*Y1));
  RCP<MV > Y2tmp1 = rcp(new MV(*Y2));
  RCP<MV > Y2tmp2 = rcp(new MV(*Y2));
  
  if (F_location_ == 1) 
  {
    // First block row
    poissonOnPoissonMatrix_->apply(*X0, *Y0);
    poissonOnDensityMatrix_->apply(*X2, *Y0tmp);
    Y0->update(1.0, *Y0tmp, 1.0);
    // Second block row
    cmsOnPoissonMatrix_->apply(*X0, *Y1);
    cmsOnCmsMatrix2_->apply(*X1, *Y1tmp1);
    cmsOnDensityMatrix_->apply(*X2, *Y1tmp2);
    Y1->update(1.0, *Y1tmp1, 1.0);
    Y1->update(1.0, *Y1tmp2, 1.0);
    // Third block row
    Y2->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *densityOnCmsMatrix_, *X1, 0.0);
    Y2->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *densityOnDensityMatrix_, *X2, 1.0);
  } //end if
  else 
  {
    // First block row
    poissonOnPoissonMatrix_->apply(*X0, *Y0);
    poissonOnDensityMatrix_->apply(*X1, *Y0tmp);
    Y0->update(1.0, *Y0tmp, 1.0);
    // Second block row
    Y1->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *densityOnDensityMatrix_, *X1, 0.0);
    Y1->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *densityOnCmsMatrix_, *X2, 1.0);
    // Third block row
    cmsOnPoissonMatrix_->apply(*X0, *Y2);
    cmsOnDensityMatrix_->apply(*X1, *Y2tmp1);
    cmsOnCmsMatrix2_->apply(*X2, *Y2tmp2);
    Y2->update(1.0, *Y2tmp1, 1.0);
    Y2->update(1.0, *Y2tmp2, 1.0);
  } //end else
} //end Apply
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
Check
(bool verbose) const 
{
  RCP<VEC > x = rcp(new VEC(getDomainMap()));
  RCP<VEC > b = rcp(new VEC(getRangeMap()));
  RCP<VEC > bb = rcp(new VEC(*b));
  x->randomize(); // Fill x with random numbers

  RCP<MV> x0 = x->offsetViewNonConst(poissonMap_, 0);
  RCP<MV> b0 = b->offsetViewNonConst(poissonMap_, 0);
  RCP<MV> x1 = x->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements());
  // Start x1 to view middle numCmsElements elements of x
  RCP<MV> b1 = b->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements());
  // Start b1 to view middle numDensity elements of b
  RCP<MV> b2 = b->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements()+densityMap_->getNodeNumElements());
  // Start b2 to view last numDensity elements of b
  RCP<MV> x2 = x->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements()+densityMap_->getNodeNumElements());
  //Start x2 to view last numCms elements of x

  // The poisson-on-poisson matrix is singular. Make x0 orthogonal to the kernel.
  RCP<MV> ones = rcp(new MV(*b0));
  //  RCP<VEC> ones = onesmv->getVectorNonConst(0);

  Scalar alpha;
  Array<Scalar> norms;
  ArrayView<Scalar> dots;
  ones->putScalar(1.0);
  ones->norm2( norms );
  alpha = norms[0];
  alpha = 1.0 / alpha;
  ones->scale( alpha );
  ones->dot( *x0, dots );
  alpha = -1.0*alpha;
  x0->update( alpha, *ones, 1.0 ); 

  apply(*x, *b); // Forward operation
  
  // Inverse if not exact, so we must modify b first:
  if (F_location_ == 1) 
  {
    b2->elementWiseMultiply(-1.0, *densityOnCmsMatrix_, *x1, 1.0);
  } //end if
  else 
  {
    b1->elementWiseMultiply(-1.0, *densityOnCmsMatrix_, *x2, 1.0);
  } //end else
  
  applyInverse(*b, *bb); // Reverse operation
  bb->update(-1.0, *x, 1.0); // Should be zero
  Scalar resid = b->norm2();
  
  if (verbose) 
  {
    std::cout << "A22 self-check residual = " << resid << std::endl;
  } //end if

  TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n"); 
} //end Check
#if LINSOLVE_PREC == 0
// Use float
template class dft_PolyA22_Coulomb_Tpetra_Operator<float, int, int>;
#elif LINSOLVE_PREC == 1
// Use double
template class dft_PolyA22_Coulomb_Tpetra_Operator<double, int, int>;
#elif LINSOLVE_PREC == 2
// Use quad double
template class dft_PolyA22_Coulomb_Tpetra_Operator<qd_real, int, int>;
#endif

