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

#ifndef TPETRA_OPERATORAPPLYINVERSE_HPP
#define TPETRA_OPERATORAPPLYINVERSE_HPP

#include "Tpetra_Operator.hpp"
#include "Tpetra_MultiVector.hpp"

namespace Tpetra{

  //! \brief A pure virtual interface for Operators with the applyInverse() method.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal, \c GlobalOrdinal and \c Node.
    The \c LocalOrdinal type, if omitted, defaults to \c int.
    The \c GlobalOrdinal type defaults to the \c LocalOrdinal type.
    The \c Node type defaults to the default node in Kokkos.
  */
  template <class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class OperatorApplyInverse : virtual public Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
  public:
    //! @name Destructor Method
    //@{

    //! Destructor.
    virtual ~OperatorApplyInverse();

    //! Apply the inverse of this operator to the multivectors.
    virtual void applyInverse(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
			      MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y) const = 0;

  }; // class OperatorApplyInverse

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  OperatorApplyInverse<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~OperatorApplyInverse()
  {
  }

} // namespace Tpetra

#endif //TPETRA_OPERATORAPPLYINVERSE_HPP
