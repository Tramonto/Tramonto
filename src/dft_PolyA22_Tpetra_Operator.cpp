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
  cmsOnCmsMatrix2_ = rcp(new MAT(cmsMap, 0));
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
  TEUCHOS_TEST_FOR_EXCEPTION(isGraphStructureSet_, std::runtime_error, "Graph structure must be set.\n");
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
  } //end if
} //end initializeProblemValues
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value, GlobalOrdinal blockColFlag )
{
  // if poisson then blockColFlag = 0
  // if density then blockColFlag = 1
  // if cms then blockColFlag = 2

  Array<GlobalOrdinal> cols(1);
  Array<Scalar> vals(1);
  cols[0] = colGID;
  vals[0] = value;

  if (cmsMap_->isNodeGlobalElement(rowGID)) { // Insert into cmsOnCmsMatrix or cmsOnDensityMatrix
    if ( blockColFlag == 2 ) { // Insert into cmsOnCmsMatrix
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
	//cout<< "row GID = " << rowGID << " value = " << value << " colGID = " << colGID << endl;
      	cmsOnDensityMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
      }
    }
    else {
      char err_msg[200];
      sprintf(err_msg,"PolyA22_Epetra_Operator::insertMatrixValue(): Invalid argument -- row in cmsMap, but blockColFlag not set for cms or density equations.");
      TEUCHOS_TEST_FOR_EXCEPT_MSG(1, err_msg);
    }
  } // end Insert into cmsOnCmsMatrix or cmsOnDensityMatrix
  else if (densityMap_->isNodeGlobalElement(rowGID)) { // Insert into densityOnDensityMatrix or densityOnCmsMatrix
    if ( blockColFlag == 1 ) { // Insert into densityOnDensityMatrix
      TEUCHOS_TEST_FOR_EXCEPT(rowGID!=colGID); // Confirm that this is a diagonal value
      densityOnDensityMatrix_->sumIntoLocalValue(densityMap_->getLocalElement(rowGID), value);
    }
    else if ( blockColFlag == 2) { // Insert into densityOnCmsMatrix
      TEUCHOS_TEST_FOR_EXCEPT(densityMap_->getLocalElement(rowGID)!=cmsMap_->getLocalElement(colGID)); // Confirm that this is a diagonal value
      densityOnCmsMatrix_->sumIntoLocalValue(densityMap_->getLocalElement(rowGID), value);
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
} //end insertMatrixValue
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  // Fill row of cmsOnCms and cmsOnDensity matrices
  if (!curRowValuesCmsOnDensity_.empty()) {
    size_t numEntriesCmsOnDensity = curRowValuesCmsOnDensity_.size();
    if (numEntriesCmsOnDensity>indicesCmsOnDensity_.size()) {
      indicesCmsOnDensity_.resize(numEntriesCmsOnDensity);
      valuesCmsOnDensity_.resize(numEntriesCmsOnDensity);
    }
    LocalOrdinal i=0;
    typename std::map<LocalOrdinal, Scalar>::iterator pos;
    for (pos = curRowValuesCmsOnDensity_.begin(); pos != curRowValuesCmsOnDensity_.end(); ++pos) {
      indicesCmsOnDensity_[i] = pos->first;
      valuesCmsOnDensity_[i++] = pos->second;
    }
    cmsOnDensityMatrix_->insertGlobalValues(curRow_, indicesCmsOnDensity_, valuesCmsOnDensity_);
  }

  if (!curRowValuesCmsOnCms_.empty()) {
    size_t numEntriesCmsOnCms = curRowValuesCmsOnCms_.size();
    if (numEntriesCmsOnCms>indicesCmsOnCms_.size()) {
      indicesCmsOnCms_.resize(numEntriesCmsOnCms);
      valuesCmsOnCms_.resize(numEntriesCmsOnCms);
    }
    LocalOrdinal i=0;
    typename std::map<LocalOrdinal, Scalar>::iterator pos;
    for (pos = curRowValuesCmsOnCms_.begin(); pos != curRowValuesCmsOnCms_.end(); ++pos) {
      indicesCmsOnCms_[i] = pos->first;
      valuesCmsOnCms_[i++] = pos->second;
    }
    cmsOnCmsMatrix2_->insertGlobalValues(curRow_, indicesCmsOnCms_, valuesCmsOnCms_);
  }

  indicesCmsOnDensity_.clear();
  valuesCmsOnDensity_.clear();
  indicesCmsOnCms_.clear();
  valuesCmsOnCms_.clear();
  curRowValuesCmsOnDensity_.clear();
  curRowValuesCmsOnCms_.clear();

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

  insertRow(); // Dump any remaining entries
  cmsOnCmsMatrix2_->fillComplete();
  if (!isFLinear_) {
    insertRow(); // Dump any remaining entries
    cmsOnDensityMatrix_->fillComplete(densityMap_, cmsMap_);
  } //end if

  if (!hasDensityOnCms)  // Confirm that densityOnCmsMatrix is zero
  {
    Scalar normvalue = densityOnCmsMatrix_->normInf();
    TEUCHOS_TEST_FOR_EXCEPT(normvalue!=0.0);
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

  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
#ifdef KDEBUG
  printf("\n\n\n\ndft_PolyA22_Tpetra_Operator::applyInverse()\n\n\n\n");
#endif

  size_t NumVectors = Y.getNumVectors();
  size_t numCmsElements = cmsMap_->getNodeNumElements();
  size_t numDensityElements = densityMap_->getNodeNumElements();

  RCP<VEC> tmp = rcp(new VEC(densityMap_));
  RCP<VEC> tmp2 = rcp(new VEC(densityMap_));

  tmp->reciprocal(*densityOnDensityMatrix_);

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

    // Second block row: Y2 = DD\X2
    Y2->elementWiseMultiply(1.0, *tmp, *X2, 0.0);
    // First block row: Y1 = CC \ (X1 - CD*Y2)
    cmsOnDensityMatrix_->apply(*Y2, *Y1tmp);
    Y1tmp->update(1.0, *X1, -1.0);
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    VEC cmsOnCmsDiag(cmsMap_);
    cmsOnCmsMatrix2_->getLocalDiagCopy(cmsOnCmsDiag);
    tmp2->reciprocal(cmsOnCmsDiag);
    Y1->elementWiseMultiply(1.0, *tmp2, *Y1tmp, 0.0);

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

    // First block row: Y1 = DD\X1
    Y1->elementWiseMultiply(1.0, *tmp, *X1, 0.0);
    // Second block row: Y2 = CC \ (X2 - CD*Y1)
    cmsOnDensityMatrix_->apply(*Y1, *Y2tmp);
    Y2tmp->update(1.0, *X2, -1.0);
    // Extract diagonal of cmsOnCmsMatrix and use that as preconditioner
    VEC cmsOnCmsDiag(cmsMap_);
    cmsOnCmsMatrix2_->getLocalDiagCopy(cmsOnCmsDiag);
    tmp2->reciprocal(cmsOnCmsDiag);
    Y2->elementWiseMultiply(1.0, *tmp2, *Y2tmp, 0.0);

  } //end else

} //end applyInverse
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{
  TEUCHOS_TEST_FOR_EXCEPT(!X.getMap()->isSameAs(*getDomainMap()));
  TEUCHOS_TEST_FOR_EXCEPT(!Y.getMap()->isSameAs(*getRangeMap()));
  TEUCHOS_TEST_FOR_EXCEPT(Y.getNumVectors()!=X.getNumVectors());
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
    RCP<MV > Y1tmp = rcp(new MV(*Y1));
    cmsOnCmsMatrix2_->apply(*X1, *Y1tmp);
    Y1->update(1.0, *Y1tmp, 1.0);
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
    RCP<MV > Y2tmp = rcp(new MV(*Y2));
    cmsOnCmsMatrix2_->apply(*X2, *Y2tmp);
    Y2->update(1.0, *Y2tmp, 1.0);

  } //end else

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
    std::cout << "A22 self-check residual = " << resid << std::endl;
  } //end if

  TEUCHOS_TEST_FOR_EXCEPTION(resid > 1.0E-12, std::runtime_error, "Bad residual.\n");

} //end Check
#if LINSOLVE_PREC == 0
// Use float
template class dft_PolyA22_Tpetra_Operator<float, int, int>;
#elif LINSOLVE_PREC == 1
// Use double
template class dft_PolyA22_Tpetra_Operator<double, int, int>;
#elif LINSOLVE_PREC == 2
// Use quad double
template class dft_PolyA22_Tpetra_Operator<qd_real, int, int>;
#elif LINSOLVE_PREC == 3
// Use double double
template class dft_PolyA22_Tpetra_Operator<dd_real, int, int>;
#endif
