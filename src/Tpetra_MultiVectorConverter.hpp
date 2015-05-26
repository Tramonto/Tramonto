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

namespace Tpetra {

  //! MultiVectorConverter: A means of converting a Tpetra Multivector to a different scalar precision.

  template<class DomainScalar, class RangeScalar=DomainScalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Details::DefaultTypes::node_type>
  class MultiVectorConverter {
  public:

    /** \name Typedefs */
    //@{
    typedef Tpetra::MultiVector<DomainScalar,LocalOrdinal,GlobalOrdinal,Node> DMV;
    typedef Tpetra::MultiVector<RangeScalar,LocalOrdinal,GlobalOrdinal,Node> RMV;

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

    inline
    void convert(const DMV& X, RMV& Y)
    {
      // Convert X from DomainScalar precision to RangeScalar precision
      for (size_t j=0; j<X.getNumVectors(); ++j) {
	ArrayRCP<const DomainScalar> xvecVals = X.getVector( j )->get1dView();
	if( xvecVals.size() ) {
	  std::transform( xvecVals.begin(), 
			  xvecVals.end(), 
			  Y.getVectorNonConst( j )->get1dViewNonConst().begin(), 
			  Teuchos::asFunc<RangeScalar>() );
	}
      }
      return;
    }

  };

}

#endif /* TPETRA_MULTIVECTORCONVERTER_H */
