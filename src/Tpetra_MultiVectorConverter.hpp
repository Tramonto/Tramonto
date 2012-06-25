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

#ifndef TPETRA_MULTIVECTORCONVERTER_HPP
#define TPETRA_MULTIVECTORCONVERTER_HPP

#include "BelosTypes.hpp"

namespace Tpetra {

  //! MultiVectorConverter: A means of converting a Tpetra Multivector to a different scalar precision.

  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MultiVectorConverter {
  public:

    /** \name Typedefs that give access to the template parameters. */
    //@{

    /// \typedef scalar_type
    /// \brief The type of the entries of the input and output multivectors.
    typedef typename Teuchos::ScalarTraits<Scalar>::halfPrecision halfScalar;
    typedef typename Teuchos::ScalarTraits<Scalar>::doublePrecision doubleScalar;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;
    typedef Tpetra::MultiVector<halfScalar,LocalOrdinal,GlobalOrdinal,Node> MV_H;
    typedef Tpetra::MultiVector<doubleScalar,LocalOrdinal,GlobalOrdinal,Node> MV_D;

    //@}

    //! @name Constructor
    //@{
    MultiVectorConverter()
    {
      return;
    }

    //! Destructor
    virtual ~MultiVectorConverter() {}
    //@}

    void scalarToHalf(const MV& X, MV_H Y)
    {
      // Demote X from scalar precision to halfPrecision
      Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > scalarVec;
      for (size_t j=0; j<X.getNumVectors(); j++) {
	scalarVec = X.getVector( j );
	ArrayRCP<const Scalar> vecVals = scalarVec->get1dView();
	Array<halfScalar> new_vals;
	new_vals.resize(vecVals.size());
	if( vecVals.size() > 0 ) {
	  std::transform( vecVals.begin(), vecVals.end(), new_vals.begin(), Teuchos::asFunc<halfScalar>() );
	}
	for( LocalOrdinal li=X.getMap()->getMinLocalIndex(); li <= X.getMap()->getMaxLocalIndex(); ++li ) {
	  Y.replaceLocalValue( li, j, new_vals[li] );
	}
      }
    }

    void doubleToScalar(const MV_D& X, MV Y)
    {
      // Demote X from doubleScalar precision to scalar precision
      Teuchos::RCP<const Tpetra::Vector<doubleScalar,LocalOrdinal,GlobalOrdinal,Node> > doubleScalarVec;
      for (size_t j=0; j<X.getNumVectors(); j++) {
	doubleScalarVec = X.getVector( j );
	ArrayRCP<const doubleScalar> vecVals = doubleScalarVec->get1dView();
	Array<Scalar> new_vals;
	new_vals.resize(vecVals.size());
	if( vecVals.size() > 0 ) {
	  std::transform( vecVals.begin(), vecVals.end(), new_vals.begin(), Teuchos::asFunc<Scalar>() );
	}
	for( LocalOrdinal li=X.getMap()->getMinLocalIndex(); li <= X.getMap()->getMaxLocalIndex(); ++li ) {
	  Y.replaceLocalValue( li, j, new_vals[li] );
	}
      }
    }

    void scalarToDouble(const MV& X, MV_D Y)
    {
      // Promote X from scalar precision to doublePrecision
      Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > scalarVec;
      for (size_t j=0; j<X.getNumVectors(); j++) {
	scalarVec = X.getVector( j );
	ArrayRCP<const Scalar> vecVals = scalarVec->get1dView();
	Array<doubleScalar> new_vals;
	new_vals.resize(vecVals.size());
	if( vecVals.size() > 0 ) {
	  std::transform( vecVals.begin(), vecVals.end(), new_vals.begin(), Teuchos::asFunc<doubleScalar>() );
	}
	for( LocalOrdinal li=X.getMap()->getMinLocalIndex(); li <= X.getMap()->getMaxLocalIndex(); ++li ) {
	  Y.replaceLocalValue( li, j, new_vals[li] );
	}
      }
    }

    void halfToScalar(const MV_H& X , MV& Y)
    {
      // Promote X from scalar precision to doublePrecision
      Teuchos::RCP<const Tpetra::Vector<halfScalar,LocalOrdinal,GlobalOrdinal,Node> > halfScalarVec;
      for (size_t j=0; j<X.getNumVectors(); j++) {
	halfScalarVec = X.getVector( j );
	ArrayRCP<const halfScalar> vecVals = halfScalarVec->get1dView();
	Array<doubleScalar> new_vals;
	new_vals.resize(vecVals.size());
	if( vecVals.size() > 0 ) {
	  std::transform( vecVals.begin(), vecVals.end(), new_vals.begin(), Teuchos::asFunc<Scalar>() );
	}
	for( LocalOrdinal li=X.getMap()->getMinLocalIndex(); li <= X.getMap()->getMaxLocalIndex(); ++li ) {
	  Y.replaceLocalValue( li, j, new_vals[li] );
	}
      }
    }

  };

}

#endif /* TPETRA_MULTIVECTORCONVERTER_H */
