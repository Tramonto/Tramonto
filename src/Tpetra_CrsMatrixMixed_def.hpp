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

#ifndef TPETRA_CRSMATRIXMIXED_DEF_HPP
#define TPETRA_CRSMATRIXMIXED_DEF_HPP

// TODO: row-wise insertion of entries in globalAssemble() may be more efficient
// TODO: consider maintaining sorted entries at all times and leaning heavily on STL set_intersect/set_union methods for all insert/replace/suminto

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_as.hpp>

#include "Tpetra_CrsMatrixMultiplyOp.hpp" // must include for implicit instantiation to work
#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_CrsMatrixMixed_decl.hpp"
#endif

namespace Tpetra {

  template <class Scalar,
	    class LocalOrdinal,
	    class GlobalOrdinal,
	    class Node,
	    class LocalMatOps>
  CrsMatrixMixed<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrixMixed (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
		  size_t maxNumEntriesPerRow,
		  ProfileType pftype,
		  const RCP<Teuchos::ParameterList>& params)
    : CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> (rowMap, maxNumEntriesPerRow, pftype, params)
  {
    typedef RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >    CRSM;
    doubleScalarMultiplyOp_ = createCrsMatrixMultiplyOp<double_scalar_type,Scalar> ((CRSM)rcp (this,false).getConst ());
  }

  template <class Scalar,
	    class LocalOrdinal,
	    class GlobalOrdinal,
	    class Node,
	    class LocalMatOps>
  CrsMatrixMixed<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrixMixed (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &rowMap,
		  const ArrayRCP<const LocalOrdinal> &NumEntriesPerRowToAlloc,
		  ProfileType pftype,
		  const RCP<Teuchos::ParameterList>& params)
    : CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> (rowMap, NumEntriesPerRowToAlloc, pftype, params)
  {
    typedef RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >    CRSM;
    doubleScalarMultiplyOp_ = createCrsMatrixMultiplyOp<double_scalar_type,Scalar> ((CRSM)rcp (this,false).getConst ());
  }

  template <class Scalar,
	    class LocalOrdinal,
	    class GlobalOrdinal,
	    class Node,
	    class LocalMatOps>
  CrsMatrixMixed<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrixMixed (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
		  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
		  size_t maxNumEntriesPerRow,
		  ProfileType pftype,
		  const RCP<Teuchos::ParameterList>& params)
    : CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> (rowMap, colMap, maxNumEntriesPerRow, pftype, params)
  {
    typedef RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >    CRSM;
    doubleScalarMultiplyOp_ = createCrsMatrixMultiplyOp<double_scalar_type,Scalar> ((CRSM)rcp (this,false).getConst ());
  }

  template <class Scalar,
	    class LocalOrdinal,
	    class GlobalOrdinal,
	    class Node,
	    class LocalMatOps>
  CrsMatrixMixed<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrixMixed (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
		  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
		  const ArrayRCP<const LocalOrdinal> &NumEntriesPerRowToAlloc,
		  ProfileType pftype,
		  const RCP<Teuchos::ParameterList>& params)
    : CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> (rowMap, colMap, NumEntriesPerRowToAlloc, pftype, params)
  {
    typedef RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >    CRSM;
    doubleScalarMultiplyOp_ = createCrsMatrixMultiplyOp<double_scalar_type,Scalar> ((CRSM)rcp (this,false).getConst ());
  }


  template<class Scalar,
	   class LocalOrdinal,
	   class GlobalOrdinal,
	   class Node,
	   class LocalMatOps>
  CrsMatrixMixed<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  CrsMatrixMixed (const RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > &graph,
		  const RCP<Teuchos::ParameterList>& params)
    : CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> (graph, params)
  {
    typedef RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >    CRSM;
    doubleScalarMultiplyOp_ = createCrsMatrixMultiplyOp<double_scalar_type,Scalar> ((CRSM)rcp (this,false).getConst ());
  }

  template<class Scalar,
	   class LocalOrdinal,
	   class GlobalOrdinal,
	   class Node,
	   class LocalMatOps>
  CrsMatrixMixed<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::
  ~CrsMatrixMixed() {}


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                  User-visible class methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void CrsMatrixMixed<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::applyDouble(
	    const MultiVector<typename Teuchos::ScalarTraits<Scalar>::doublePrecision,LocalOrdinal,GlobalOrdinal,Node> &X,
	    MultiVector<typename Teuchos::ScalarTraits<Scalar>::doublePrecision,LocalOrdinal,GlobalOrdinal,Node> &Y,
	    Teuchos::ETransp mode,
	    typename Teuchos::ScalarTraits<Scalar>::doublePrecision alpha,
	    typename Teuchos::ScalarTraits<Scalar>::doublePrecision beta) const {
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() == false, std::runtime_error,
				typeName(*this) << "::apply(): cannot call apply() until fillComplete() has been called.");
    doubleScalarMultiplyOp_->apply(X,Y,mode,alpha,beta);
  }

} // namespace Tpetra

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_CRSMATRIXMIXED_INSTANT(SCALAR,LO,GO,NODE) \
  \
  template class CrsMatrixMixed< SCALAR , LO , GO , NODE >;

#define TPETRA_CRSMATRIXMIXED_CONVERT_INSTANT(SO,SN,LO,GO,NODE) \
  \
  template RCP< CrsMatrixMixed< SN , LO , GO , NODE > >   \
		CrsMatrixMixed< SO , LO , GO , NODE >::convert< SN > () const;

#define TPETRA_CRSMATRIXMIXED_IMPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  RCP<CrsMatrixMixed<SCALAR, LO, GO, NODE> >                                \
  importAndFillCompleteCrsMatrixMixed (const RCP<const CrsMatrixMixed<SCALAR, LO, GO, NODE> >& sourceMatrix, \
				  const Import<CrsMatrixMixed<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
					       CrsMatrixMixed<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
					       CrsMatrixMixed<SCALAR, LO, GO, NODE>::node_type>& importer, \
				  const RCP<const Map<CrsMatrixMixed<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
				  const RCP<const Map<CrsMatrixMixed<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
							       const RCP<Teuchos::ParameterList>& params);

#define TPETRA_CRSMATRIXMIXED_EXPORT_AND_FILL_COMPLETE_INSTANT(SCALAR, LO, GO, NODE) \
  template<>                                                                        \
  RCP<CrsMatrixMixed<SCALAR, LO, GO, NODE> >                                \
  exportAndFillCompleteCrsMatrixMixed (const RCP<const CrsMatrixMixed<SCALAR, LO, GO, NODE> >& sourceMatrix, \
				  const Export<CrsMatrixMixed<SCALAR, LO, GO, NODE>::local_ordinal_type,  \
					       CrsMatrixMixed<SCALAR, LO, GO, NODE>::global_ordinal_type,  \
					       CrsMatrixMixed<SCALAR, LO, GO, NODE>::node_type>& exporter, \
				  const RCP<const Map<CrsMatrixMixed<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::node_type> >& domainMap, \
				  const RCP<const Map<CrsMatrixMixed<SCALAR, LO, GO, NODE>::local_ordinal_type,      \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::global_ordinal_type,     \
							       CrsMatrixMixed<SCALAR, LO, GO, NODE>::node_type> >& rangeMap,  \
							       const RCP<Teuchos::ParameterList>& params);

#endif
