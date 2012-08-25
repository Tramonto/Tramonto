
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

#ifndef DFT_A22MATRIX_TPETRA_OPERATOR_H
#define DFT_A22MATRIX_TPETRA_OPERATOR_H

#include "Tpetra_Headers.hpp"

//! dft_A22Matrix_Epetra_Operator: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! Special 2*numBeads by 2*numBeads for Tramonto polymer problems.
*/

template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node=Kokkos::DefaultNode::DefaultNodeType>
class dft_A22Matrix_Tpetra_Operator:
  public virtual Tpetra::OperatorApplyInverse<Scalar,LocalOrdinal,GlobalOrdinal,Node>
{

public:
  TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node);
  TYPEDEF_MIXED(Scalar, LocalOrdinal, GlobalOrdinal, Node);

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads system

  dft_A22Matrix_Tpetra_Operator
  (const RCP<const MAP > & block2Map, RCP<ParameterList> parameterList);
  //@}
  //@{ \name Assembly methods.
  void
  initializeProblemValues
  ();

  void
  insertMatrixValue
  (GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value);

  void
  finalizeProblemValues
  ();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_A22Matrix_Tpetra_Operator();
  //@}

  //@{ \name Atribute set methods.

    //! Unsupported feature, returns -1.
  void
  SetUseTranspose
  (bool UseTranspose)
  {
    return;
  };
  //@}

  //@{ \name Mathematical functions.

    //! Returns the result of a dft_A22Matrix_Epetra_Operator applied to a Epetra_MultiVector X in Y.
    /*!
    \param In
	   X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  virtual void
  apply
  (const MV& X, MV& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, Scalar alpha = 1.0, Scalar beta = 0.0) const;

  //! Returns the result of an inverse dft_A22Matrix_Epetra_Operator applied to a Epetra_MultiVector X in Y.
  /*!
    \param In
    X - An Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y - An Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
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
  NormInf() const
  {
    return(0.0);
  };
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
  UseTranspose() const
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

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const RCP<const COMM> &
  Comm
  () const
  {
    return(block2Map_->getComm());
  };

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const RCP<const MAP> &
  getDomainMap
  () const
  {
    return(block2Map_);
  };

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const RCP<const MAP> &
  getRangeMap
  () const
  {
    return(block2Map_);
  };

  //! Returns a pointer to the Epetra_CrsMatrix object that is the A22 matrix
  RCP<MAT_P>
  getA22Matrix
  ()
  {
    return(A22Matrix_);
  }
  //@}

private:

  void insertRow();
  const RCP<const MAP > block2Map_;
  RCP<ParameterList> parameterList_;
  RCP<MAT_P> A22Matrix_;
  RCP<PRECOND_P> A22Inverse_;
  RCP<MOP> A22InverseMixed_;
  const char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool firstTime_;
  GlobalOrdinal curRow_;
  std::map<GlobalOrdinal, precScalar> curRowValues_;
  Array<GlobalOrdinal> indices_;
  Array<precScalar> values_;
  Ifpack2::Factory factory_;
};

#endif /* DFT_A22MATRIX_TPETRA_OPERATOR_H */
