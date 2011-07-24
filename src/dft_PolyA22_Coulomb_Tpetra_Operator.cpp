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
  
  poissonMatrix_ = rcp(new MAT(poissonMap, 0));
  cmsOnPoissonMatrix_ = rcp(new MAT(cmsMap, 0));
  poissonOnDensityMatrix_ = rcp(new MAT(poissonMap, 0));
  Label_ = "dft_PolyA22_Coulomb_Tpetra_Operator";
  cmsOnDensityMatrix_->setObjectLabel("PolyA22Coulomb::cmsOnDensityMatrix");
  poissonMatrix_->setObjectLabel("PolyA22Coulomb::poissonMatrix");
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
      cmsOnDensityMatrix_->setAllToScalar(0.0); //not needed
    } //end if
    cmsOnCmsMatrix_->putScalar(0.0);
    densityOnDensityMatrix_->putScalar(0.0);
    densityOnCmsMatrix_->putScalar(0.0);
    poissonMatrix_->setAllToScalar(0.0); //not needed
    cmsOnPoissonMatrix_->setAllToScalar(0.0); //not needed
    poissonOnDensityMatrix_->setAllToScalar(0.0); //not needed
  } //end if
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value) 
{

  Array<GlobalOrdinal> cols(1);
  cols[0] = colGID;
  Array<Scalar> vals(1);
  vals[0] = value;

  /* The poissonMatrix_, poissonOnDensityMatrix_, cmsOnPoissonMatrix_, and cmsOnDensityMatrix_ values do not change between iterations */  

  if (poissonMap_->getLocalElement(rowGID)) 
  {
    if (poissonMap_->getLocalElement(colGID)) 
    {
      if (firstTime_) 
      { 
	if (rowGID != curPoissonRow_) 
        {
	  insertPoissonRow();
	  curPoissonRow_ = rowGID;
	} //end if
	curPoissonRowValues_[colGID] += value;
      }
      else 
      {
	poissonMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
      } //end else
    } //end if
    else  //!poissonMap_.MyGID(colGID)
    {
      if (firstTime_) 
      {
	if (rowGID != curPDRow_) 
        {
	  insertPDRow();
	  curPDRow_ = rowGID;
	} //end if
	curPDRowValues_[colGID] += value;
      } //end if
      else 
      {
	poissonOnDensityMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
      } //end else
    } //end else
  } //end if
  else if (cmsMap_->getLocalElement(rowGID)) 
  {
    if (rowGID==colGID)
    {
      cmsOnCmsMatrix_->sumIntoGlobalValue(rowGID, value); // Storing this cms block in a vector since it is diagonal
    } //end if
    else if (densityMap_->getLocalElement(colGID))
    {
      if (firstTime_) 
      {
	if (rowGID!=curRow_) 
        { 
	  P22TO::insertRow();
          // Dump the current contents of curRowValues_ into matrix and clear map
	  curRow_=rowGID;
	} //end if
	curRowValues_[colGID] += value;
      } //end if
      else if (!isFLinear_) 
      {
	GlobalOrdinal newRowGID = densityMap_->getGlobalElement(cmsMap_->getLocalElement(rowGID));
	cmsOnDensityMatrix_->sumIntoGlobalValues(newRowGID, cols, vals);
      } //end elseif
    } //end elseif
    else if (poissonMap_->getLocalElement(colGID)) 
    {
      if (firstTime_) 
      {
	if (rowGID!=curCPRow_) 
        {
	  insertCPRow(); 
          // Dump the current contents of curCPRowValues_ into matrix and clear map
	  curCPRow_=rowGID;
	} //end if
	curCPRowValues_[colGID] += value;
      } //end if
      else 
      {
	cmsOnPoissonMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
      } //end else
    } //end elseif
  } //end if
  else 
  {
    if (rowGID==colGID)
      densityOnDensityMatrix_->sumIntoGlobalValue(rowGID, value); // Storing this density block in a vector since it is diagonal
    else
    {
      densityOnCmsMatrix_->sumIntoGlobalValue(rowGID, value); 
      // Storing this density block in a vector since it is diagonal
    } //end else
  } //end else
} //end insertMatrixValue
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertPoissonRow
() 
{
  if (curPoissonRowValues_.empty()) 
  {
    return;
  } //end if
  size_t numEntries = curPoissonRowValues_.size();
  if (numEntries>indices_.size()) 
  {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  } //end if
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos=curPoissonRowValues_.begin(); pos!=curPoissonRowValues_.end(); ++pos) 
  {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  } //end for
  poissonMatrix_->insertGlobalValues(curPoissonRow_, indices_, values_);
  curPoissonRowValues_.clear();
} //end insertPoissonRow
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertCPRow
() 
{
  if (curCPRowValues_.empty()) 
  {
    return;
  } //end if
  size_t numEntries = curCPRowValues_.size();
  if (numEntries>indices_.size()) 
  {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  } //end if
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos=curCPRowValues_.begin(); pos!=curCPRowValues_.end(); ++pos) 
  {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  } //end for
  cmsOnPoissonMatrix_->insertGlobalValues(curCPRow_, indices_, values_);
  curCPRowValues_.clear();
} //end insertCPRow
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertPDRow
() 
{
  if (curPDRowValues_.empty()) 
  {
    return;
  } //end if
  size_t numEntries = curPDRowValues_.size();
  if (numEntries>indices_.size()) 
  {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  } //end if
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos=curPDRowValues_.begin(); pos!=curPDRowValues_.end(); ++pos) 
  {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  } //end for
  poissonOnDensityMatrix_->insertGlobalValues(curPDRow_, indices_, values_); 
  curPDRowValues_.clear();
} //end insertPDRow
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
  
  if (firstTime_) 
  {
    insertPoissonRow();
    insertCPRow();
    insertPDRow();
    if (!isFLinear_) 
    {
      P22TO::insertRow(); 
      // Dump any remaining entries
    } //end if
  } //end if
  
  if (!isFLinear_) 
  {
    cmsOnDensityMatrix_->fillComplete();
  } //end if
  
  //only need to do the following the firstTime_
  poissonMatrix_->fillComplete();
  cmsOnPoissonMatrix_->fillComplete(poissonMap_, cmsMap_);
  poissonOnDensityMatrix_->fillComplete(densityMap_, poissonMap_);

  /* This appears to make minimal difference, so could be removed */
  for (LocalOrdinal i = 0; i < poissonMap_->getNodeNumElements(); i++) 
  {
    GlobalOrdinal row = poissonMatrix_->getRowMap()->getGlobalElement(i);
    Array<GlobalOrdinal> indices(1);
    indices[0] = row;
    Array<Scalar> values(1);
    values[0] = 10.0e-12;
    poissonMatrix_->sumIntoGlobalValues(row, indices, values);
  } //end for
  
  if (firstTime_) 
  {
    ////////////////MLPrec = rcp(new ML_Epetra::MultiLevelPreconditioner(*poissonMatrix_, *MLList_, true));
  } //end if
  
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
  RCP<MV > Y1tmp = rcp(new MV(*Y1));
  RCP<MV > Y1tmp2 = rcp(new MV(*Y1));
  RCP<MV > Y2tmp = rcp(new MV(*Y2));
  RCP<MV > Y2tmp2 = rcp(new MV(*Y2));
  
  RCP<VEC> tmp;
  densityOnDensityMatrix_->reciprocal(*tmp);
  if (F_location_ == 1) 
  {
    Y2->elementWiseMultiply(1.0, *tmp, *X2, 0.0);
    poissonOnDensityMatrix_->apply(*Y2, *Y0tmp);
    Y0->update(1.0, *X0, -1.0, *Y0tmp, 0.0);
  } //end if
  else 
  {
    Y1->elementWiseMultiply(1.0, *tmp, *X1, 0.0);
    poissonOnDensityMatrix_->apply(*Y1, *Y0tmp);
    Y0->update(1.0, *X0, -1.0, *Y0tmp, 0.0);
  } //end else
    
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

  if (F_location_ == 1) 
  {
    cmsOnDensityMatrix_->apply(*Y2, *Y1tmp);
    cmsOnPoissonMatrix_->apply(*Y0, *Y1tmp2);
    Y1tmp->update(1.0, *Y1tmp2, 1.0);
    Y1->update(1.0, *X1, -1.0, *Y1tmp, 0.0);
    RCP<VEC> tmp2;
    cmsOnCmsMatrix_->reciprocal(*tmp2);
    Y1->elementWiseMultiply(1.0, *tmp2, *Y1, 0.0);
  } //end if
  else 
  { 
    cmsOnDensityMatrix_->apply(*Y1, *Y2tmp);
    cmsOnPoissonMatrix_->apply(*Y0, *Y2tmp2);
    Y2tmp->update(1.0, *Y2tmp2, 1.0);
    Y2->update(1.0, *X2, -1.0, *Y2tmp, 0.0);
    RCP<VEC> tmp;
    cmsOnCmsMatrix_->reciprocal(*tmp);
    Y2->elementWiseMultiply(1.0, *tmp, *Y2, 0.0);
  } //end else
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
  RCP<MV > Y1tmp = rcp(new MV(*Y1));
  RCP<MV > Y1tmp2 = rcp(new MV(*Y1));
  RCP<MV > Y2tmp = rcp(new MV(*Y2));
  RCP<MV > Y2tmp2 = rcp(new MV(*Y2));
  
  poissonMatrix_->apply(*X0, *Y0tmp);
  
  if (F_location_ == 1) 
  {
    poissonOnDensityMatrix_->apply(*X2, *Y0tmp2);
    cmsOnPoissonMatrix_->apply(*X0, *Y1tmp);
    cmsOnDensityMatrix_->apply(*X2, *Y1tmp2);
    Y2tmp->elementWiseMultiply(1.0, *densityOnCmsMatrix_, *X1, 0.0);
    Y2tmp2->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X2, 0.0);
    Y1->elementWiseMultiply(1.0, *cmsOnCmsMatrix_, *X1, 0.0);
    Y1->update(1.0, *Y1tmp, 1.0,  *Y1tmp2, 1.0);
    Y2->update(1.0, *Y2tmp, 1.0, *Y2tmp2, 0.0);
    Y0->update(1.0, *Y0tmp, 1.0, *Y0tmp2, 0.0);
  } //end if
  else 
  {
    poissonOnDensityMatrix_->apply(*X1, *Y0tmp2);
    cmsOnPoissonMatrix_->apply(*X0, *Y2tmp);
    cmsOnDensityMatrix_->apply(*X1, *Y2tmp2);
    Y1tmp->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X1, 0.0);
    Y1tmp2->elementWiseMultiply(1.0, *densityOnCmsMatrix_, *X2, 0.0);
    Y2->elementWiseMultiply(1.0, *cmsOnCmsMatrix_, *X2, 0.0);
    Y2->update(1.0, *Y2tmp, 1.0, *Y2tmp2, 1.0);
    Y1->update(1.0, *Y1tmp, 1.0, *Y1tmp2, 0.0);
    Y0->update(1.0, *Y0tmp, 1.0, *Y0tmp2, 0.0);
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
  x->randomize(); // Fill x with random numbers

  RCP<MV> x1 = x->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements());
  // Start x1 to view middle numCmsElements elements of x
  RCP<MV> b2 = b->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements()+cmsMap_->getNodeNumElements());
  // Start b2 to view last numDensity elements of b
  RCP<MV> x2 = x->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements()+densityMap_->getNodeNumElements());
  //Start x2 to view last numCms elements of x
  RCP<MV> b1 = b->offsetViewNonConst(densityMap_, poissonMap_->getNodeNumElements());
  // Start b1 to view middle numDensity elements of b

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
  
  applyInverse(*b, *b); // Reverse operation
  b->update(-1.0, *x, 1.0); // Should be zero
  Scalar resid = b->norm2();
  
  if (verbose) 
  {
    std::cout << "A22 self-check residual = " << resid << endl;
  } //end if

  TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n"); 
} //end Check
template class dft_PolyA22_Coulomb_Tpetra_Operator<double, int, int>;
