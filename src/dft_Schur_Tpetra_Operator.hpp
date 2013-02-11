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

#ifndef DFT_SCHUR_TPETRA_OPERATOR_H
#define DFT_SCHUR_TPETRA_OPERATOR_H

#include "Tpetra_Headers.hpp"
#include "dft_PolyA11_Tpetra_Operator.hpp"
#include "dft_PolyA11_Coulomb_Tpetra_Operator.hpp"
#include "dft_PolyA22_Tpetra_Operator.hpp"
#include "dft_PolyA22_Coulomb_Tpetra_Operator.hpp"

/*! Special 2-by-2 block operator for Tramonto polymer and explicit non-local density problems.
*/

template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
  class Node=Kokkos::DefaultNode::DefaultNodeType>
class dft_Schur_Tpetra_Operator:
  public virtual Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{

 public:
TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node)
TYPEDEF_MIXED(Scalar, LocalOrdinal, GlobalOrdinal, Node)

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2-by-2 block system

  dft_Schur_Tpetra_Operator
  (RCP<APINV> A11, RCP<MAT_P> A12, RCP<MAT_P> A21, RCP<APINV> A22);
  //@}
  //@{ \name Destructor.
    //! Destructor
  ~dft_Schur_Tpetra_Operator();
  //@}

  //@{ \name Atribute set methods.
  void
  SetSchurComponents
  (RCP<MAT_P> A11invMatrix, RCP<MAT_P> A22Matrix)
  {
    A11invMatrix_ = A11invMatrix;
    A22Matrix_ = A22Matrix;
  }

    //! Unsupported feature, throws an exception.
  void
  SetUseTranspose
  (bool UseTranspose)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "SetUseTranspose is not supported.\n");
  };
  //@}

  //@{ \name Mathematical functions.

    //! Returns the result of a dft_Schur_Tpetra_Operator applied to a MultiVector X in Y.
    /*!
    \param In
	   X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A MultiVector of dimension NumVectors containing result.
  */
  void
  apply
  (const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const;

  //! Unsupported feature, throws an exception.
  void
  applyInverse
  (const MV& X, MV& Y) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "applyInverse is not supported.\n");
  };

  //! Generate RHS.
  void
  ComputeRHS
  (const MV& B1, const MV& B2, MV& B) const;

  //! Generate X1.
  void
  ComputeX1
  (const MV& B1, const MV& X2, MV& X1) const;

  //! Apply global operator.
  void
  ApplyGlobal
  (const MV& X1, const MV& X2, MV& Y1, MV& Y2) const;

  //! Explicitly form Schur complement as an CrsMatrix and return a pointer to it.
  void formSchurComplement();

  //! Returns a pointer to the Tpetra_CrsMatrix object that is the schur complement
  RCP<MAT_P>
  getSchurComplement
  ()
  {
    if (A11invMatrix_.is_null() || A22Matrix_.is_null())
      {
	RCP<MAT_P> null;
	return(null);  // We cannot form S without the component matrices
      } //end if
    formSchurComplement();
    return(S_);
  }

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
  const RCP<const COMM> &
  Comm
  () const
  {
    return(A22_->getDomainMap()->getComm());
  };

  //! Returns the Map object associated with the domain of this operator.
  const RCP<const MAP > &
  getDomainMap
  () const
  {
    return(A22_->getDomainMap());
  };

  //! Returns the Map object associated with the range of this operator.
  const RCP<const MAP>&
  getRangeMap
  () const
  {
    return(A22_->getRangeMap());
  };
  //@}

  RCP<APINV> A11_;
  /*!< The 1,1 block of the 2 by 2 block matrix */
  RCP<MAT_P> A12_;
  RCP<MMOP_P> A12op_;
  /*!< The 1,2 block of the 2 by 2 block matrix */
  RCP<MAT_P> A21_;
  RCP<MMOP_P> A21op_;
  /*!< The 2,1 block of the 2 by 2 block matrix */
  RCP<APINV> A22_;
  /*!< The 2,2 block of the 2 by 2 block matrix */
  RCP<MAT_P> A11invMatrix_;
  /*!< The inverse of A11 in matrix form, if available */
  RCP<MAT_P> A22Matrix_;
  /*!< A22 as a matrix, if available */
  RCP<MAT_P> A11invA12_;
  /* Intermediate matrix containing inv(A11)*A12 */
  RCP<MAT_P> A21A11invA12_;
  /* Intermediate matrix containing A21*inv(A11)*A12 */
  RCP<MAT_P> S_;
  /* Schur complement, if formed */
  const char * Label_; /*!< Description of object */
};

#endif /* DFT_SCHUR_TPETRA_OPERATOR_H */
