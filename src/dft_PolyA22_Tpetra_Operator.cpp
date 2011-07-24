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

#include "dft_PolyA22_Tpetra_Operator.hpp"

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA22_Tpetra_Operator
(const RCP<const MAP > & cmsMap, const RCP<const MAP > & densityMap, const RCP<const MAP > & block2Map, 
 RCP<ParameterList> parameterList) 
  : cmsMap_(cmsMap),
    densityMap_(densityMap),
    block2Map_(block2Map),
    parameterList_(parameterList),
    Label_(0),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    isFLinear_(false),
    firstTime_(true),
    curRow_(-1) 
{

  cmsOnDensityMatrix_ = rcp(new MAT(cmsMap, 0));
  cmsOnCmsMatrix_ = rcp(new VEC(cmsMap));
  densityOnDensityMatrix_ = rcp(new VEC(densityMap));
  densityOnCmsMatrix_ = rcp(new VEC(densityMap));
  Label_ = "dft_PolyA22_Tpetra_Operator";
  cmsOnDensityMatrix_->setObjectLabel("PolyA22::cmsOnDensityMatrix");
  F_location_ = Teuchos::getParameter<LocalOrdinal>(*parameterList_, "F_location"); 
  //F in NE if F_location = 1, F in SW otherwise
} //end constructor
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_PolyA22_Tpetra_Operator
() 
{
} //end destructor
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
() 
{
  TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n"); 
  isLinearProblemSet_ = false; // We are reinitializing the linear problem

  if (!firstTime_) 
  {
    if (!isFLinear_) 
    {
      cmsOnDensityMatrix_->setAllToScalar(0.0);
    } //end if
    cmsOnCmsMatrix_->putScalar(0.0);
    densityOnDensityMatrix_->putScalar(0.0);
    densityOnCmsMatrix_->putScalar(0.0);
  } //end if
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value) 
{

  Array<GlobalOrdinal> cols(1);
  Array<Scalar> vals(1);
  cols[0] = colGID;
  vals[0] = value;
  if (cmsMap_->isNodeGlobalElement(rowGID)) 
  {
    if (rowGID==colGID)
    {
      cmsOnCmsMatrix_->sumIntoGlobalValue(rowGID, value); 
      // Storing this cms block in a vector since it is diagonal
    } //end if
    else
    {
      if (firstTime_) 
      {
	if (rowGID!=curRow_) 
        {
	  insertRow();  
          // Dump the current contents of curRowValues_ into matrix and clear map
	  curRow_=rowGID;
	} //end if
	curRowValues_[colGID] += value;
      } //end if
      else if (!isFLinear_) 
      { 
	cmsOnDensityMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
      } //end elseif
    } //end else
  } //end if
  else 
  {
    if (rowGID==colGID)
    {
      densityOnDensityMatrix_->sumIntoGlobalValue(rowGID, value); 
      // Storing this density block in a vector since it is diagonal
    } //end if
    else 
    {
      TEST_FOR_EXCEPT(densityMap_->getLocalElement(rowGID)!=cmsMap_->getLocalElement(colGID)); 
      // Confirm that this is a diagonal value
      densityOnCmsMatrix_->sumIntoGlobalValue(rowGID, value); 
      // Storing this density block in a vector since it is diagonal
    } //end else
  } //end else
} //end insertMatrixValue
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
() 
{
  if (curRowValues_.empty())
  {
    return;
  } //end if
  size_t numEntries = curRowValues_.size();
  if (numEntries>indices_.size()) 
  {
    indices_.resize(numEntries);
    values_.resize(numEntries);
  } //end if
  LocalOrdinal i=0;
  typename std::map<GlobalOrdinal, Scalar>::iterator pos;
  for (pos = curRowValues_.begin(); pos != curRowValues_.end(); ++pos) 
  {
    indices_[i] = pos->first;
    values_[i++] = pos->second;
  } //end for
  
  cmsOnDensityMatrix_->insertGlobalValues(curRow_, indices_, values_);

  curRowValues_.clear();
} //end insertRow
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
() 
{
  if (isLinearProblemSet_) 
  {
    return; // nothing to do
  } //end if
  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_->getGlobalNumElements()==densityMap_->getGlobalNumElements(); 

  if (!isFLinear_) 
  {
    insertRow(); // Dump any remaining entries
    cmsOnDensityMatrix_->fillComplete(densityMap_, cmsMap_);
  } //end if

  if (!hasDensityOnCms)  // Confirm that densityOnCmsMatrix is zero
  {
    Scalar normvalue = densityOnCmsMatrix_->normInf();
    TEST_FOR_EXCEPT(normvalue!=0.0);
  } //end if
  //cout << "CmsOnDensityMatrix Inf Norm = " << cmsOnDensityMatrix_->NormInf() << endl;
  //densityOnDensityMatrix_->NormInf(&normvalue);
  //cout << "DensityOnDensityMatrix Inf Norm = " << normvalue << endl;
  //cout << "DensityOnCmsMatrix Inf Norm = " << normvalue << endl;

  isLinearProblemSet_ = true;
  firstTime_ = false;
} //end finalizeProblemValues
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const 
{
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
  
  // replacing Ddc with a zero matrix for the applyInverse method only.

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
  
  // replacing Ddc with a zero matrix for the applyInverse method only.

  // Our algorithm is then:
  // Y1 = Ddd \ X1
  // Y2 = Dcc \ (X2 - F*Y1)  

  //Scalar normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in PolyA22 applyInverse = " << normvalue << endl;

  TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap())); 
  TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  printf("\n\n\n\ndft_PolyA22_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  size_t NumVectors = Y.getNumVectors();
  size_t numCmsElements = cmsMap_->getNodeNumElements();
  size_t numDensityElements = densityMap_->getNodeNumElements();

  if (F_location_ == 1)  //F in NE
  {
    RCP<MV > Y1 = Y.offsetViewNonConst(cmsMap_, 0);
    // Y1 is a view of the first numCms elements of Y
    RCP<MV > Y2 = Y.offsetViewNonConst(densityMap_, numCmsElements);
    // Start Y2 to view last numDensity elements of Y
    RCP<const MV > X1 = X.offsetView(cmsMap_, 0);
    // Start X1 to view first numCms elements of X
    RCP<const MV > X2 = X.offsetView(densityMap_, numCmsElements);
    // Start X2 to view last numDensity elements of X
    RCP<MV > Y1tmp = rcp(new MV(*Y1));
    RCP<VEC> temp1;
    densityOnDensityMatrix_->reciprocal(*temp1);
    Y2->elementWiseMultiply(1.0, *temp1, *X2, 0.0);
    cmsOnDensityMatrix_->apply(*Y2, *Y1tmp);
    Y1->update(1.0, *X1, -1.0, *Y1tmp, 0.0);
    RCP<VEC> temp2;
    cmsOnCmsMatrix_->reciprocal(*temp2);
    Y1->elementWiseMultiply(1.0, *temp2, *Y1, 0.0);
  } //end if
  else 
  {
    RCP<MV > Y1 = Y.offsetViewNonConst(densityMap_, 0);
    // Y1 is a view of the first numDensity elements of Y
    RCP<MV > Y2 = Y.offsetViewNonConst(cmsMap_, numDensityElements);
    // Start Y2 to view last numCms elements of Y
    RCP<const MV > X1 = X.offsetView(densityMap_, 0);
    // Start X1 to view first numDensity elements of X
    RCP<const MV > X2 = X.offsetView(cmsMap_, numDensityElements);
    // Start X2 to view last numCms elements of X
    RCP<MV > Y2tmp = rcp(new MV(*Y2));
    RCP<VEC> temp1;
    densityOnDensityMatrix_->reciprocal(*temp1);
    Y1->elementWiseMultiply(1.0, *temp1, *X1, 0.0);
    cmsOnDensityMatrix_->apply(*Y1, *Y2tmp);
    Y2->update(1.0, *X2, -1.0, *Y2tmp, 0.0);
    RCP<VEC> temp2;
    cmsOnCmsMatrix_->reciprocal(*temp2);
    Y2->elementWiseMultiply(1.0, *temp2, *Y2, 0.0);
  } //end else
  
  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA22 applyInverse = " << normvalue << endl;
} //end applyInverse
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const 
{
  //Scalar normvalue;
  //X.NormInf(&normvalue);
  //cout << "Norm of X in PolyA22 Apply = " << normvalue << endl;

  TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
  size_t NumVectors = Y.getNumVectors();
  size_t numCmsElements = cmsMap_->getNodeNumElements();
  size_t numDensityElements = densityMap_->getNodeNumElements(); 

  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_->getGlobalNumElements()==densityMap_->getGlobalNumElements(); 

  if (F_location_ == 1) //F in NE
  {
    RCP<MV> Y1 = Y.offsetViewNonConst(cmsMap_, 0);
    // Y1 is a view of the first numCms elements of Y
    RCP<MV> Y2 = Y.offsetViewNonConst(densityMap_, numCmsElements);
    // Start Y2 to view last numDensity elements of Y
    RCP<const MV> X1 = X.offsetView(cmsMap_, 0);
    // Start X1 to view first numCms elements of X
    RCP<const MV> X2 = X.offsetView(densityMap_, numCmsElements);
    // Start X2 to view last numDensity elements of X

    cmsOnDensityMatrix_->apply(*X2, *Y1);
    Y1->elementWiseMultiply(1.0, *cmsOnCmsMatrix_, *X1, 1.0);
    Y2->elementWiseMultiply(1.0, *densityOnCmsMatrix_, *X1, 0.0);
    Y2->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X2, 1.0);
  } //end if
  else
  {
  
    RCP<MV> Y1 = Y.offsetViewNonConst(densityMap_, 0);//rcp(new MV(densityMap_, Y1ptr, NumVectors));
    // Y1 is a view of the first numDensity elements of Y
    RCP<MV> Y2 = Y.offsetViewNonConst(cmsMap_, numDensityElements);//rcp(new MV(cmsMap_, Y2ptr, NumVectors)); 
    // Start Y2 to view last numCms elements of Y
    RCP<const MV> X1 = X.offsetView(densityMap_, 0);//rcp(new MV(densityMap_, X1ptr, NumVectors));
    // Start X1 to view first numDensity elements of X
    RCP<const MV> X2 = X.offsetView(cmsMap_, numDensityElements);//rcp(new MV(cmsMap_, X2ptr, NumVectors)); 
    // Start X2 to view last numCms elements of X
    if (hasDensityOnCms) 
    {
      // convert X2 map
      RCP<const MV> X2tmp = X.offsetView(densityMap_, numDensityElements);
      Y1->elementWiseMultiply(1.0, *densityOnCmsMatrix_, *X2tmp, 0.0); 
      // Momentarily make X2 compatible with densityMap_
    } //end if
    Y1->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X1, 1.0);
    cmsOnDensityMatrix_->apply(*X1, *Y2);
    Y2->elementWiseMultiply(1.0, *cmsOnCmsMatrix_, *X2, 1.0);
  } //end else
  

  //Y.NormInf(&normvalue);
  //cout << "Norm of Y in PolyA22 Apply = " << normvalue << endl;
} //end Apply
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
Check
(bool verbose) const 
{
  RCP<VEC > x = rcp(new VEC(getDomainMap()));
  RCP<VEC > b = rcp(new VEC(getRangeMap()));
  x->randomize(); // Fill x with random numbers

  // densityOnCmsMatrix will be nonzero only if cms and density maps are the same size
  bool hasDensityOnCms = cmsMap_->getGlobalNumElements()==densityMap_->getGlobalNumElements(); 

  apply(*x, *b); // Forward operation

  if (hasDensityOnCms) 
  {
    if (F_location_ == 1) //F in NE
    {
      // Inverse is not exact, so we must modify b2 first:
      RCP<MV> x1 = x->offsetViewNonConst(densityMap_, 0);
      // Start x1 to view first numCmsElements elements of x
      RCP<MV> b2 = b->offsetViewNonConst(densityMap_, cmsMap_->getNodeNumElements());
      // Start b2 to view last numDensity elements of b
      b2->elementWiseMultiply(-1.0, *densityOnCmsMatrix_, *x1, 1.0);
    } //end if
    else 
    {
      // Inverse is not exact, so we must modify b1 first:
      RCP<MV> x2 = x->offsetViewNonConst(densityMap_, densityMap_->getNodeNumElements());
      //Start x2 to view last numCms elements of x
      RCP<MV> b1 = b->offsetViewNonConst(densityMap_, 0);
      // Start b1 to view first numDensity elements of b
      b1->elementWiseMultiply(-1.0, *densityOnCmsMatrix_, *x2, 1.0);
    } //end else
  } //end if

  applyInverse(*b, *b); // Reverse operation

  b->update(-1.0, *x, 1.0); // Should be zero

  Scalar resid = b->norm2();

  if (verbose) 
  {
    std::cout << "A22 self-check residual = " << resid << endl;
  } //end if

  TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n");
} //end Check
template class dft_PolyA22_Tpetra_Operator<double, int, int>;
