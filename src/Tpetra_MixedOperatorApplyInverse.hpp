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

#ifndef TPETRA_MIXEDOPERATORAPPLYINVERSE_HPP
#define TPETRA_MIXEDOPERATORAPPLYINVERSE_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_OperatorApplyInverse.hpp"
#include "Tpetra_MultiVectorConverter.hpp"

#define SCALAR_P      0
#define HALF_SCALAR_P 1

namespace Tpetra {

  //! MixedOperatorApplyInverse: An implementation of the OperatorApplyInverse class that does the Operator apply() and applyInverse() to a doublePrecision vector.
  /*! The MixedOperatorApplyInverse class implements OperatorApplyInverse using another pre-constructed OperatorApplyInverse object.
    Once constructed, an MixedOperatorApplyInverse can be used to apply the input operator to doublePrecision vectors
    as long as the appropriate apply and applyInverse methods are implemented in the original OperatorApplyInverse object.
  */
  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MixedOperatorApplyInverse: public virtual Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
  public:

    /** \name Typedefs that give access to the template parameters. */
    //@{
    typedef typename Teuchos::ScalarTraits<Scalar>::halfPrecision halfScalar;

    //@}

    //! @name Constructor
    //@{
    //! Uses an Operator instance to implement the Operator interface.
    /*! Facilitates the use of an Operator instance as an inverse operator.
      \param In - A fully-constructed Operator object.
    */
    MixedOperatorApplyInverse(Teuchos::RCP<Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node> > operatorIn) {
      operator_ = operatorIn;
      operatorPrec_ = SCALAR_P;
      return;
    }
#if MIXED_PREC == 1
    MixedOperatorApplyInverse(Teuchos::RCP<Tpetra::OperatorApplyInverse<halfScalar, LocalOrdinal, GlobalOrdinal, Node> > operatorIn) {
      operatorHalf_ = operatorIn;
      operatorPrec_ = HALF_SCALAR_P;
      return;
    }
#endif

    //! Destructor
    virtual ~MixedOperatorApplyInverse(){}
    //@}

    //! @name Mathematical functions
    //@{

    //! Returns the result of an MixedOperatorApplyInverse applied to a MultiVector X in Y.
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
      if (operatorPrec_ == HALF_SCALAR_P) {

	Teuchos::RCP<Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node> > halfX =
	  Teuchos::rcp( new Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node>(X.getMap(), X.getNumVectors()) );
	Teuchos::RCP<Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node> > halfY =
	  Teuchos::rcp( new Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node>(Y.getMap(), Y.getNumVectors()) );
	// Cast alpha, beta to halfPrecision
	halfScalar alphaHalf = Teuchos::as<halfScalar>(alpha);
	halfScalar betaHalf = Teuchos::as<halfScalar>(beta);

	// Demote X from scalar precision to halfPrecision
	mvConverter->scalarToHalf( X, *halfX );

	// Apply the operator
	operatorHalf_->apply(*halfX,*halfY,mode,alphaHalf,betaHalf);

	// Promote halfY from halfPrecision to scalar precision
	mvConverter->halfToScalar( *halfY, Y );

      } else if (operatorPrec_ == SCALAR_P) {

	// Apply the operator
	operator_->apply(X,Y,mode,alpha,beta);

      }

    }

    //! Returns the result of an MixedOperatorApplyInverse inverse applied to a MultiVector X in Y.

    void applyInverse(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
		      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
    {
      if (operatorPrec_ == HALF_SCALAR_P) {

	Teuchos::RCP<Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node> > halfX =
	  Teuchos::rcp( new Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node>( X.getMap(), X.getNumVectors() ) );
	Teuchos::RCP<Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node> > halfY =
	  Teuchos::rcp( new Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node>( Y.getMap(), Y.getNumVectors() ) );

	// Demote X from scalar precision to halfPrecision
	mvConverter->scalarToHalf( X, *halfX );

	// Apply the inverse operator
	operatorHalf_->applyInverse(*halfX,*halfY);

	// Promote halfY from halfPrecision to scalar precision
	mvConverter->halfToScalar( *halfY, Y );

      } else if (operatorPrec_ == SCALAR_P) {

	// Apply the operator
	operator_->applyInverse(X,Y);

      }

    }

  //! @name Attribute access functions
  //@{

  //! Returns a pointer to the Operator operator object that was used to create this MixedOperatorApplyInverse object.
    Teuchos::RCP<Tpetra::OperatorApplyInverse<halfScalar, LocalOrdinal, GlobalOrdinal, Node> > getOperator() const {return(operator_);}

  //! Returns the BlockMap object associated with the domain of this matrix operator.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const
    {
      if (operatorPrec_ == HALF_SCALAR_P) {
	return(operatorHalf_->getDomainMap());
      } else if (operatorPrec_ == SCALAR_P) {
	return(operator_->getDomainMap());
      }
    }

  //! Returns the BlockMap object associated with the range of this matrix operator.
    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const
    {
      if (operatorPrec_ == HALF_SCALAR_P) {
	return(operatorHalf_->getRangeMap());
      } else if (operatorPrec_ == SCALAR_P) {
	return(operator_->getRangeMap());
      }
    }
  //@}

 protected:
    int operatorPrec_;
    Teuchos::RCP<Tpetra::OperatorApplyInverse<Scalar, LocalOrdinal, GlobalOrdinal, Node> > operator_;
    Teuchos::RCP<Tpetra::OperatorApplyInverse<halfScalar, LocalOrdinal, GlobalOrdinal, Node> > operatorHalf_;
    Teuchos::RCP<Tpetra::MultiVectorConverter<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mvConverter;
};

}

#endif /* TPETRA_MIXEDOPERATORAPPLYINVERSE_H */
