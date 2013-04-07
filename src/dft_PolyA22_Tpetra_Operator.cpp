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

  cmsOnDensityMatrix_ = rcp(new MAT_P(cmsMap, 0));
  cmsOnDensityMatrixOp_ = rcp(new MMOP_P(cmsOnDensityMatrix_));
  cmsOnCmsMatrix_ = rcp(new MAT_P(cmsMap, 0));
  cmsOnCmsMatrixOp_ = rcp(new MMOP_P(cmsOnCmsMatrix_));
  densityOnDensityMatrix_ = rcp(new VEC(densityMap));
  densityOnDensityInverse_ = rcp(new VEC(densityMap));
  densityOnCmsMatrix_ = rcp(new MAT_P(densityMap, 0));
  densityOnCmsMatrixOp_ = rcp(new MMOP_P(densityOnCmsMatrix_));
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
    cmsOnCmsMatrix_->resumeFill();
    cmsOnCmsMatrix_->setAllToScalar(0.0);
    densityOnDensityMatrix_->putScalar(0.0);
    densityOnDensityInverse_->putScalar(0.0);
    densityOnCmsMatrix_->resumeFill();
    densityOnCmsMatrix_->setAllToScalar(0.0);
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
  Array<precScalar> vals(1);
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
	cmsOnCmsMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
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
      // The density-on-cms matrix is diagonal for most but not all use cases.
      if (firstTime_) {
	if (rowGID!=curRow_) {
	  insertRow();  // Dump the current contents of curRowValues maps  into matrix and clear map
	  curRow_=rowGID;
	}
	curRowValuesDensityOnCms_[colGID] += value;
      }
      else
      	densityOnCmsMatrix_->sumIntoGlobalValues(rowGID, cols, vals);
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
    ITER pos;
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
    ITER pos;
    for (pos = curRowValuesCmsOnCms_.begin(); pos != curRowValuesCmsOnCms_.end(); ++pos) {
      indicesCmsOnCms_[i] = pos->first;
      valuesCmsOnCms_[i++] = pos->second;
    }
    cmsOnCmsMatrix_->insertGlobalValues(curRow_, indicesCmsOnCms_, valuesCmsOnCms_);
  }

  if (!curRowValuesDensityOnCms_.empty()) {
    size_t numEntriesDensityOnCms = curRowValuesDensityOnCms_.size();
    if (numEntriesDensityOnCms>indicesDensityOnCms_.size()) {
      indicesDensityOnCms_.resize(numEntriesDensityOnCms);
      valuesDensityOnCms_.resize(numEntriesDensityOnCms);
    }
    LocalOrdinal i=0;
    ITER pos;
    for (pos = curRowValuesDensityOnCms_.begin(); pos != curRowValuesDensityOnCms_.end(); ++pos) {
      indicesDensityOnCms_[i] = pos->first;
      valuesDensityOnCms_[i++] = pos->second;
    }
    densityOnCmsMatrix_->insertGlobalValues(curRow_, indicesDensityOnCms_, valuesDensityOnCms_);
  }

  indicesCmsOnDensity_.clear();
  valuesCmsOnDensity_.clear();
  indicesCmsOnCms_.clear();
  valuesCmsOnCms_.clear();
  indicesDensityOnCms_.clear();
  valuesDensityOnCms_.clear();
  curRowValuesCmsOnDensity_.clear();
  curRowValuesCmsOnCms_.clear();
  curRowValuesDensityOnCms_.clear();

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
  RCP<ParameterList> pl = rcp(new ParameterList(parameterList_->sublist("fillCompleteList")));
  cmsOnCmsMatrix_->fillComplete(pl);
  if (!isFLinear_) {
    insertRow(); // Dump any remaining entries
    cmsOnDensityMatrix_->fillComplete(densityMap_, cmsMap_, pl);
  } //end if

  if (!hasDensityOnCms)  // Confirm that densityOnCmsMatrix is zero
  {
    //    Scalar normvalue = densityOnCmsMatrix_->normInf();
    //    TEUCHOS_TEST_FOR_EXCEPT(normvalue!=0.0);
  } else {
    insertRow(); // Dump any remaining entries
    densityOnCmsMatrix_->fillComplete(cmsMap_, densityMap_, pl);
  }

  // Form the inverse of the densityOnDensityMatrix
  densityOnDensityInverse_->reciprocal(*densityOnDensityMatrix_);

  // Use a diagonal preconditioner for the cmsOnCmsMatrix
  RCP<const MAT_P> const_matrix = Teuchos::rcp_implicit_cast<const MAT_P>(cmsOnCmsMatrix_);
  cmsOnCmsInverse_ = rcp(new PRECOND_D(const_matrix));
  cmsOnCmsInverseOp_ = rcp(new PRECOND_D_OP(cmsOnCmsInverse_));
  TEUCHOS_TEST_FOR_EXCEPT(cmsOnCmsInverse_==Teuchos::null);
  cmsOnCmsInverse_->initialize();
  cmsOnCmsInverse_->compute();

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

  if (F_location_ == 1)
  {
    //F in NE

    // Y1 is a view of the first numCms elements of Y
    RCP<MV > Y1 = Y.offsetViewNonConst(cmsMap_, 0);
    // Y2 is a view of the last numDensity elements of Y
    RCP<MV > Y2 = Y.offsetViewNonConst(densityMap_, numCmsElements);
    // X1 is a view of the first numCms elements of X
    RCP<const MV > X1 = X.offsetView(cmsMap_, 0);
    // X2 is a view of the last numDensity elements of X
    RCP<const MV > X2 = X.offsetView(densityMap_, numCmsElements);
    RCP<MV > Y1tmp = rcp(new MV(*Y1));

    // Second block row: Y2 = DD\X2
    Y2->elementWiseMultiply(1.0, *densityOnDensityInverse_, *X2, 0.0);

    // First block row: Y1 = CC \ (X1 - CD*Y2)
    cmsOnDensityMatrixOp_->apply(*Y2, *Y1tmp);
    Y1tmp->update(1.0, *X1, -1.0);
    cmsOnCmsInverseOp_->apply(*Y1tmp, *Y1);

  }
  else
  {
    //F in SW

    // Y1 is a view of the first numDensity elements of Y
    RCP<MV > Y1 = Y.offsetViewNonConst(densityMap_, 0);
    // Y2 is a view of the last numCms elements of Y
    RCP<MV > Y2 = Y.offsetViewNonConst(cmsMap_, numDensityElements);
    // X1 is a view of the first numDensity elements of X
    RCP<const MV > X1 = X.offsetView(densityMap_, 0);
    // X2 is a view of the last numCms elements of X
    RCP<const MV > X2 = X.offsetView(cmsMap_, numDensityElements);
    RCP<MV > Y2tmp = rcp(new MV(*Y2));

    // First block row: Y1 = DD\X1
    Y1->elementWiseMultiply(1.0, *densityOnDensityInverse_, *X1, 0.0);

    // Second block row: Y2 = CC \ (X2 - CD*Y1)
    cmsOnDensityMatrixOp_->apply(*Y1, *Y2tmp);
    Y2tmp->update(1.0, *X2, -1.0);
    cmsOnCmsInverseOp_->apply(*Y2tmp, *Y2);
  }

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

  if (F_location_ == 1)
  {
    //F in NE

    // Y1 is a view of the first numCms elements of Y
    RCP<MV> Y1 = Y.offsetViewNonConst(cmsMap_, 0);
    // Y2 is a view of the last numDensity elements of Y
    RCP<MV> Y2 = Y.offsetViewNonConst(densityMap_, numCmsElements);
    // X1 is a view of the first numCms elements of X
    RCP<const MV> X1 = X.offsetView(cmsMap_, 0);
    // X2 is a view of the last numDensity elements of X
    RCP<const MV> X2 = X.offsetView(densityMap_, numCmsElements);

    // First block row
    cmsOnDensityMatrixOp_->apply(*X2, *Y1);
    RCP<MV > Y1tmp = rcp(new MV(*Y1));
    cmsOnCmsMatrixOp_->apply(*X1, *Y1tmp);
    Y1->update(1.0, *Y1tmp, 1.0);

    // Second block row
    if (hasDensityOnCms) {
      densityOnCmsMatrixOp_->apply(*X1, *Y2);
      Y2->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X2, 1.0);
    } else {
      Y2->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X2, 0.0);
    }

  }
  else
  {
    //F in SW

    // Y1 is a view of the first numDensity elements of Y
    RCP<MV> Y1 = Y.offsetViewNonConst(densityMap_, 0);
    // Y2 is a view of the last numCms elements of Y
    RCP<MV> Y2 = Y.offsetViewNonConst(cmsMap_, numDensityElements);
    // X1 is a view of the first numDensity elements of X
    RCP<const MV> X1 = X.offsetView(densityMap_, 0);
    // X2 is a view of the last numCms elements of X
    RCP<const MV> X2 = X.offsetView(cmsMap_, numDensityElements);

    // First block row
    if (hasDensityOnCms) {
      densityOnCmsMatrixOp_->apply(*X2, *Y1);
      Y1->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X1, 1.0);
    } else {
      Y1->elementWiseMultiply(1.0, *densityOnDensityMatrix_, *X1, 0.0);
    }

    // Second block row
    cmsOnDensityMatrixOp_->apply(*X1, *Y2);
    RCP<MV > Y2tmp = rcp(new MV(*Y2));
    cmsOnCmsMatrixOp_->apply(*X2, *Y2tmp);
    Y2->update(1.0, *Y2tmp, 1.0);

  }

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
      RCP<MV> x1 = x->offsetViewNonConst(cmsMap_, 0);
      // Start x1 to view first numCmsElements elements of x
      RCP<MV> b2 = b->offsetViewNonConst(densityMap_, cmsMap_->getNodeNumElements());
      // Start b2 to view last numDensity elements of b
      RCP<MV > DCx1 = rcp(new MV(*b2));
      densityOnCmsMatrixOp_->apply(*x1, *DCx1);
      b2->update(-1.0, *DCx1, 1.0); // b2 = b2 - DC*x1
    }
    else
    {
      // Inverse is not exact, so we must modify b1 first:
      RCP<MV> x2 = x->offsetViewNonConst(cmsMap_, densityMap_->getNodeNumElements());
      //Start x2 to view last numCms elements of x
      RCP<MV> b1 = b->offsetViewNonConst(densityMap_, 0);
      // Start b1 to view first numDensity elements of b
      RCP<MV > DCx2 = rcp(new MV(*b1));
      densityOnCmsMatrixOp_->apply(*x2, *DCx2);
      b1->update(-1.0, *DCx2, 1.0); // b1 = b1 - DC*x2
    }
  }

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
