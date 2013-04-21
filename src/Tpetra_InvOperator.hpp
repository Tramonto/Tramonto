// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef TPETRA_INVOPERATOR_HPP
#define TPETRA_INVOPERATOR_HPP

#include "Tpetra_OperatorApplyInverse.hpp"

namespace Tpetra {

  //! InvOperator: An implementation of the OperatorApplyInverse class that reverses the role of apply() and applyInverse() methods.
  /*! The InvOperator class implements OperatorApplyInverse using another pre-constructed OperatorApplyInverse object.
    Once constructed, an InvOperator can be used as the inverse of the input operator
    object as long as the appropriate apply and applyInverse methods are implemented in the original OperatorApplyInverse object.
  */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class InvOperator: public virtual Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:

    /** \name Typedefs that give access to the template parameters. */
    //@{

    /// \typedef scalar_type
    /// \brief The type of the entries of the input and output multivectors.
    typedef Scalar scalar_type;

    /// \typedef local_ordinal_type
    /// \brief The local index type.
    typedef LocalOrdinal local_ordinal_type;

    /// \typedef global_ordinal_type
    /// \brief The global index type.
    typedef GlobalOrdinal global_ordinal_type;

    /// \typedef node_type
    /// \brief The Kokkos Node type.
    typedef Node node_type;

    //@}

    //! @name Constructor
    //@{
    //! Uses an OperatorApplyInverse instance to implement the OperatorApplyInverse interface.
    /*! Facilitates the use of an OperatorApplyInverse instance as an inverse operator.
      \param In - A fully-constructed OperatorApplyInverse object.
    */
    InvOperator(Teuchos::RCP<Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node> > operatorIn) {
      operator_ = operatorIn;
      return;
    }

    //! Destructor
    virtual ~InvOperator(){}
    //@}

    //! @name Mathematical functions
    //@{

    //! Returns the result of an InvOperator applied to a MultiVector X in Y.
    /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
      vary according to the values of \c alpha and \c beta. Specifically
      - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
      - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
    */
    void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
	       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
	       Teuchos::ETransp mode = Teuchos::NO_TRANS,
	       Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
	       Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const
    {
      operator_->applyInverse(X,Y);
    }

    //! Returns the result of an InvOperator inverse applied to a MultiVector X in Y.

    void applyInverse(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
    {
      operator_->apply(X,Y);
    }

  //! @name Attribute access functions
  //@{

  //! Returns a pointer to the OperatorApplyInverse operator object that was used to create this InvOperator object.
    Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() const {return(operator_);}

  //! Returns the BlockMap object associated with the domain of this matrix operator.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const
    {
      return(operator_->getDomainMap());
    }

  //! Returns the BlockMap object associated with the range of this matrix operator.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const
    {
      return(operator_->getRangeMap());
    }
  //@}

 protected:
    Teuchos::RCP<Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node> > operator_;
};

}

#endif /* TPETRA_INVOPERATOR_H */
