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

#ifndef DFT_HARDSPHEREA11_TPETRA_OPERATOR_H
#define DFT_HARDSPHEREA11_TPETRA_OPERATOR_H

#include "Tpetra_Headers.hpp"

//! dft_HardSphereA11_Epetra_Operator: An implementation of the Epetra_Operator class for Tramonto Schur complements.
/*! Special 2*numBeads by 2*numBeads for Tramonto polymer problems.
*/

template<class Scalar,class MatScalar=Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node = Kokkos::DefaultNode::DefaultNodeType>
class dft_HardSphereA11_Tpetra_Operator:
  public virtual Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{

public:
  TYPEDEF(Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node)

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads system

  dft_HardSphereA11_Tpetra_Operator<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>
  (const RCP<const MAP> & indNonLocalRowMap, const RCP<const MAP> & depNonLocalRowMap, const RCP<const MAP> & block1Map,
   RCP<ParameterList> parameterList);
  //@}
  //@{ \name Assembly methods.
  virtual void
  initializeProblemValues
  ();

  virtual void
  insertMatrixValue
  (LocalOrdinal rowGID, LocalOrdinal colGID, MatScalar value);

  virtual void
  finalizeProblemValues
  ();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_HardSphereA11_Tpetra_Operator();
  //@}

  //@{ \name Atribute set methods.

    //! Unsupported feature, returns -1.
  void
  SetUseTranspose
  (bool UseTranspose)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "SetUseTranspose is unsupported.\n");
  };
  //@}

  //@{ \name Mathematical functions.

    //! Returns the result of a dft_HardSphereA11_Epetra_Operator applied to a Epetra_MultiVector X in Y.
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

  //! Returns the result of an inverse dft_HardSphereA11_Epetra_Operator applied to a Epetra_MultiVector X in Y.
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
  NormInf
  () const {
    return(0.0);
  };

  //! Check for inconsistencies in operators.
  /* \param verbose (In) Print the residual of inv(A11)*A11*x_random.

     \return Returns 0 if residual is "small", otherwise it returns -1.
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

  //! Return the number of entries in this block
  virtual GlobalOrdinal
  getNumEntries
  () const
  {
    return(nnz_);
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const RCP<const COMM> &
  Comm
  () const
  {
    return(block1Map_->getComm());
  };

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual RCP<const MAP>
  getDomainMap
  () const
  {
    return(block1Map_);
  };

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual RCP<const MAP>
  getRangeMap
  () const
  {
    return(block1Map_);
  };

  //! Returns a pointer to the Epetra_CrsMatrix object that is the A22 matrix
  RCP<MAT>
  getA11invMatrix
  ()
  {
    formA11invMatrix();
    return(A11invMatrix_);
  }
  //@}

protected:

  void insertRow();
  const RCP<const MAP> indNonLocalMap_;
  const RCP<const MAP> depNonLocalMap_;
  const RCP<const MAP> block1Map_;
  RCP<MAT> matrix_;
  RCP<MMOP> matrixOperator_;
  RCP<MAT> A11invMatrix_;
  RCP<MMOP> A11invMatrixOperator_;
  GlobalOrdinal nnz_;
  const char * Label_; /*!< Description of object */
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool firstTime_;
  GlobalOrdinal curRow_;
  std::map<GlobalOrdinal, MatScalar> curRowValues_;
  Array<GlobalOrdinal> indices_;
  Array<MatScalar> values_;
  void formA11invMatrix();
  RCP<ParameterList> parameterList_;
};

#endif /* DFT_HARDSPHEREA11_TPETRA_OPERATOR_H */
