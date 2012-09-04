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

#ifndef DFT_POLYA22_TPETRA_OPERATOR_H
#define DFT_POLYA22_TPETRA_OPERATOR_H

#include "Tpetra_Headers.hpp"
#include "Ifpack.h"

//! dft_PolyA22_Tpetra_Operator: An implementation of the Tpetra_Operator class for Tramonto Schur complements.
/*! Special 2*numBeads by 2*numBeads for Tramonto polymer problems.
*/

template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
class dft_PolyA22_Tpetra_Operator:
  public virtual Tpetra::OperatorApplyInverse<Scalar,LocalOrdinal,GlobalOrdinal,Node>
{

 public:
TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node);
TYPEDEF_MIXED(Scalar, LocalOrdinal, GlobalOrdinal, Node);

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads system
  /* dft_PolyA22_Tpetra_Operator(const Map & cmsMap, const Map & densityMap, const Map & block2Map, LocalOrdinal * options, Scalar * params);*/

  dft_PolyA22_Tpetra_Operator
  (const RCP<const MAP > & cmsMap, const RCP<const MAP > & densityMap,
   const RCP<const MAP > & block2Map, RCP<ParameterList> parameterList);
  //@}
  //@{ \name Assembly methods.

  //! Assert that field dependence on primitive densities is linear; manager will not reset values between nonlinear solves.
  /*! This method can be called to assert that the field variable dependence on primitive densities does change from one linear solve to the next.
      In this case, we can avoid filling the associated matrix coefficients.  Calling this method with "true" will cause the problem manager not reset
      the matrix coefficients for this block and to ignore any values that are submitted for entry in this block.
     \param isLinear (In) Set to true if the field dependence is linear on primitive densities.
     \warning This method can be called at any time, but should be called before the initializeValues() method is called for the second solve; By default the manager assumes that the relationship is non-linear, so values will be reset to zero and must be refilled before each linear solve.
  */
  void
  setFieldOnDensityIsLinear
  (bool isLinear)
  {
    isFLinear_ = isLinear;
  }

  virtual void
  initializeProblemValues
  ();

  virtual void
  insertMatrixValue
  (GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value, GlobalOrdinal blockColFlag);

  virtual void
  finalizeProblemValues
  ();

  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_PolyA22_Tpetra_Operator();
  //@}

  //@{ \name Atribute get methods.

  //! Returns an Operator pointer that is actually the \e this object, since this class implements Operator.
  virtual RCP<OP>
  getA22Inv
  ()
  {
    return(rcp(this));
  }
  //@}

  //@{ \name Atribute set methods.

    //! Unsupported feature, throws an exception.
  void
  SetUseTranspose
  (bool UseTranspose){
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "SetUseTranspose is not supported.\n");
  };
  //@}

  //@{ \name Mathematical functions.

    //! Returns the result of a dft_PolyA22_Tpetra_Operator applied to a MultiVector X in Y.
    /*!
    \param In
	   X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A MultiVector of dimension NumVectors containing result.
  */
  virtual void
  apply
  (const MV& X, MV& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, Scalar alpha = 1.0, Scalar beta = 0.0) const;

  //! Returns the result of an inverse dft_PolyA22_Tpetra_Operator applied to a MultiVector X in Y.
  /*!
    \param In
    X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y - A MultiVector of dimension NumVectors containing result.
  */
  virtual void
  applyInverse
  (const MV& X, MV& Y) const;


  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].

     \warning This method must not be called unless HasNormInf() returns true.
  */
  Scalar
  NormInf
  () const
  {
    return(0.0);
  };

  //! Check for inconsistencies in operators.
  /* \param verbose (In) Print the residual of inv(A22)*A22*x_random.

  //! Throws an exception if the residual error is "large".
  */
  virtual void
  Check
  (bool verbose) const;

  //@}

  //@{ \name Atribute access functions

  //! Returns a character string describing the operator
  const char *
  Label
  () const
  {
    return(Label_);
  };

  //! Returns the current UseTranspose setting.
  bool
  UseTranspose
  () const
  {
    return(false);
  };

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool
  HasNormInf
  () const
  {
    return(false);
  };

  //! Returns a pointer to the Comm communicator associated with this operator.
  virtual const RCP<const COMM> &
  Comm
  () const
  {
    return(block2Map_->getComm());
  };

  //! Returns the Map object associated with the domain of this operator.
  virtual const RCP<const MAP > &
  getDomainMap
  () const
  {
    return(block2Map_);
  };

  //! Returns the Map object associated with the range of this operator.
  virtual const RCP<const MAP > &
  getRangeMap
  () const
  {
    return(block2Map_);
  };
  //@}

protected:

  LocalOrdinal F_location_;
  void insertRow();
  const RCP<const MAP > cmsMap_;
  const RCP<const MAP > densityMap_;
  const RCP<const MAP > block2Map_;
  RCP<ParameterList> parameterList_;
  RCP<MAT_P> cmsOnDensityMatrix_;
  RCP<DMOP_P> cmsOnDensityMatrixOp_;
  RCP<MAT_P> cmsOnCmsMatrix_;
  RCP<DMOP_P> cmsOnCmsMatrixOp_;
  RCP<VEC > densityOnDensityMatrix_;
  RCP<VEC > densityOnCmsMatrix_;
  char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool isFLinear_;
  bool firstTime_;
  GlobalOrdinal curRow_;
  std::map<GlobalOrdinal, precScalar> curRowValuesCmsOnDensity_, curRowValuesCmsOnCms_;
  Array<GlobalOrdinal> indicesCmsOnDensity_, indicesCmsOnCms_;
  Array<precScalar> valuesCmsOnDensity_, valuesCmsOnCms_;
  Teuchos::ParameterList IFList_;
  RCP<PRECOND> IFPrec;
  string IFPrecType; // incomplete LU
  GlobalOrdinal IFOverlapLevel;
  std::map<GlobalOrdinal, precScalar> curRowValues_;
  Array<GlobalOrdinal> indices_;
  Array<precScalar> values_;
}; //class dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>
#endif /* DFT_POLYA22_TPETRA_OPERATOR_H */
