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

#ifndef TPETRA_MULTIVECTORRECIPROCALMULTIPLY_HPP
#define TPETRA_MULTIVECTORRECIPROCALMULTIPLY_HPP

#include "BelosTypes.hpp"
#include <functional>
#include <algorithm>

namespace Tpetra {

  //! MultiVectorReciprocalMultiply: A means of element-wise multiplying a MultiVector by the reciprocal of a vector.

  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class MultiVectorReciprocalMultiply {
  public:

    /** \name Typedefs that give access to the template parameters. */
    //@{
    typedef Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> VEC;
    typedef Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> MV;

    //@}

    //! @name Constructor
    //@{
    MultiVectorReciprocalMultiply(Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > mvIn)
    {
      mv_ = mvIn;
      return;
    }

    //! Destructor
    virtual ~MultiVectorReciprocalMultiply() {}
    //@}


    void reciprocalMultiply(const VEC& V, const MV& X)
    {
      TEUCHOS_TEST_FOR_EXCEPT(mv_->getNumVectors()!=X.getNumVectors());
      TEUCHOS_TEST_FOR_EXCEPT(mv_->getMap()->getNodeNumElements()!=V.getMap()->getNodeNumElements());
      TEUCHOS_TEST_FOR_EXCEPT(V.getMap()->getNodeNumElements()!=X.getMap()->getNodeNumElements());

      Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mvVec;
      Teuchos::RCP<const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xVec;
      for (size_t j=0; j<mv_->getNumVectors(); j++) {
	mvVec = mv_->getVectorNonConst( j );
	xVec = X.getVector( j );
	ArrayRCP<Scalar> mvVecVals = mvVec->get1dViewNonConst();
	ArrayRCP<const Scalar> vVals = V.get1dView();
	ArrayRCP<const Scalar> xVecVals = xVec->get1dView();

	if ( mvVecVals.size() > 0 ) {
	  std::transform( xVecVals.begin(), xVecVals.end(), vVals.begin(), mvVecVals.begin(), std::divides<Scalar>() );
	}
      }
    }

  protected:
    Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > mv_;

};



}

#endif /* TPETRA_MULTIVECTORRECIPROCALMULTIPLY_H */
