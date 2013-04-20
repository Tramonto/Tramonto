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

#ifndef DFT_POLYA11_TPETRA_OPERATOR_H
#define DFT_POLYA11_TPETRA_OPERATOR_H

#include "Tpetra_Headers.hpp"

//! dft_PolyA11_Tpetra_Operator: An implementation of the Operator class for Tramonto Schur complements.
/*! Special 2*numBeads by 2*numBeads for Tramonto polymer problems.
*/
template<class Scalar,class MatScalar=Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node = Kokkos::DefaultNode::DefaultNodeType>
class dft_PolyA11_Tpetra_Operator:
  public virtual Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{

 public:
  TYPEDEF(Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node)

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads system

  dft_PolyA11_Tpetra_Operator<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>(const RCP<const MAP> & ownedMap, const RCP<const MAP> & block1Map, RCP<ParameterList> parameterList);
  //@}
  //@{ \name Assembly methods.
  virtual void
  initializeProblemValues
  ();

  virtual void
  insertMatrixValue
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, GlobalOrdinal rowGID, GlobalOrdinal colGID, MatScalar value);

  virtual void
  finalizeProblemValues
  ();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_PolyA11_Tpetra_Operator();
  //@}

  //@{ \name Atribute set methods.

    //! Unsupported feature, throws an exception.
  void
  SetUseTranspose
  (bool UseTranspose)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "SetUseTranspose is unsupported.\n");
  };
  //@}

  //@{ \name Mathematical functions.

    //! Returns the result of a dft_PolyA11_Tpetra_Operator applied to a MultiVector X in Y.
    /*!
    \param In
	   X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A MultiVector of dimension NumVectors containing result.
  */
  virtual void
  apply
  (const MV& X, MV& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, Scalar alpha = 1.0, Scalar beta = 0.0) const;

  //! Returns the result of an inverse dft_PolyA11_Tpetra_Operator applied to a MultiVector X in Y.
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
  /* \param verbose (In) Print the residual of inv(A11)*A11*x_random.

  //! Throws an exception if the residual error is "large".
  */
  void
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
    return(block1Map_->getComm());
  };

  //! Returns the Map object associated with the domain of this operator.
  virtual const RCP<const MAP> &
  getDomainMap
  () const
  {
    return(block1Map_);
  };

  //! Returns the Map object associated with the range of this operator.
  virtual const RCP<const MAP> &
  getRangeMap
  () const
  {
    return(block1Map_);
  };
  //@}

protected:

  void insertRow();
  const RCP<const MAP> ownedMap_;
  const RCP<const MAP> block1Map_;
  size_t numBlocks_;
  Array<RCP<MAT> > matrix_;
  Array<RCP<MMOP> > matrixOperator_;
  RCP<VEC> diagonal_;
  RCP<VEC> invDiagonal_;
  const char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  GlobalOrdinal curRow_;
  LocalOrdinal curOwnedPhysicsID_;
  LocalOrdinal curOwnedNode_;
  std::map<GlobalOrdinal, MatScalar> curRowValues_;
  Array<GlobalOrdinal> indices_;
  Array<MatScalar> values_;
  bool firstTime_;
  // RN_20100326: This following is needed for Krylov solvers in the
  // applyInverse method.
  RCP<ParameterList> parameterList_;
};
#endif /* DFT_POLYA11_TPETRA_OPERATOR_H */
