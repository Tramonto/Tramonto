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

#ifndef TPETRA_SCALINGCRSMATRIX_HPP
#define TPETRA_SCALINGCRSMATRIX_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Teuchos_Describable.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

namespace Tpetra {

  //! InvOperator: A class to scale a CrsMatrix.

  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ScalingCrsMatrix {
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
    //!
    /*!
      \param In - A fully-constructed CrsMatrix object.
    */
    ScalingCrsMatrix(Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > matrixIn) {
      matrix_ = matrixIn;
      return;
    }

    //! Destructor
    virtual ~ScalingCrsMatrix(){}
    //@}

    //! @name Mathematical functions
    //@{

    //! Returns the row scaling factors for the matrix

    LocalOrdinal getRowScaleFactors( RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > x, GlobalOrdinal norm ) const
    {
      // Matrix must be filled first
      LocalOrdinal ierr = 0;

      // Zero out target vector
      x->putScalar( Teuchos::ScalarTraits<Scalar>::zero() );
      ArrayRCP<Scalar> xp = x->getDataNonConst(0);

      if( (matrix_->getDomainMap())->isSameAs(*(x->getMap())) ) {
	for( LocalOrdinal i = Teuchos::OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(i) < matrix_->getNodeNumRows(); i++ ) {
	  LocalOrdinal NumEntries = matrix_->getNumEntriesInLocalRow( i );
	  ArrayView<const LocalOrdinal> rowIndices;
	  ArrayView<const Scalar> rowValues;
	  Scalar scale = Teuchos::ScalarTraits<Scalar>::zero();
	  Scalar rval = Teuchos::ScalarTraits<Scalar>::zero();
	  matrix_->getLocalRowView( i, rowIndices, rowValues );

	  // Compute row scale factors R(i)
	  if( Teuchos::as<size_t>(norm) == Teuchos::OrdinalTraits<GlobalOrdinal>::zero() ) {
	    // R(i) is infinity norm of row i
	    for( LocalOrdinal j = Teuchos::OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(j) < NumEntries; j++ ) {
	      rval = Teuchos::ScalarTraits<Scalar>::magnitude( rowValues[j] );
	      if( rval > scale )
		scale = rval;
	    }
	  } else if( Teuchos::as<size_t>(norm) == Teuchos::ScalarTraits<GlobalOrdinal>::one() ) {
	    // R(i) is 1-norm of row i
	    for( LocalOrdinal j = Teuchos::OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(j) < NumEntries; j++ )
	      scale += Teuchos::ScalarTraits<Scalar>::magnitude( rowValues[j] );
	  } else {
	    // R(i) is 2-norm of row i
	    for( LocalOrdinal j = Teuchos::OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(j) < NumEntries; j++ )
	      scale += Teuchos::ScalarTraits<Scalar>::pow( rowValues[j], 2 );
	  }
	  rowIndices = Teuchos::null;
	  rowValues = Teuchos::null;

	  // Invert scale factors
	  if( scale < Teuchos::ScalarTraits<Scalar>::rmin() ) {
	    if( scale == Teuchos::ScalarTraits<Scalar>::zero() )
	      ierr = 1; // Set error to 1 to signal that zero rowsum found (supercedes ierr = 2)
	    else if ( ierr != 1 )
	      ierr = 2;
	    xp[i] = Teuchos::ScalarTraits<Scalar>::rmax();
	  } else {
	    xp[i] = Teuchos::ScalarTraits<Scalar>::pow( scale, -1 );
	  }
	}
      } else if( matrix_->getRangeMap()->isSameAs(*(x->getMap())) ) {
	return(-2); // Don't handle this case for now
      } else { // x.Map different than both matrix_->DomainMap() and matrix_->RangeMap()
	return(-2);
      }

      return(0);
    }

    //! Returns the result of a CrsMatrix scaled by a vector.

    LocalOrdinal leftScale( const RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > x ) const
    {
      // Matrix must be filled
      LocalOrdinal ierr = 0;
      ArrayRCP<const Scalar> xp = x->getData(0);

      if( (matrix_->getDomainMap())->isSameAs(*(x->getMap())) ) {
	for( LocalOrdinal i = Teuchos::OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(i) < matrix_->getNodeNumRows(); i++ ) {

	  LocalOrdinal  numEntries = matrix_->getNumEntriesInLocalRow( i );
	  ArrayView<const LocalOrdinal> rowIndices;
	  ArrayView<const Scalar> rowValues;
	  Array<Scalar> scaledRowValues;
	  matrix_->getLocalRowView( i, rowIndices, rowValues );
	  Scalar scaleValue = xp[i];
	  scaledRowValues.resize(numEntries);
	  for( LocalOrdinal j = Teuchos::OrdinalTraits<LocalOrdinal>::zero(); Teuchos::as<size_t>(j) < numEntries; j++ ) {
	    scaledRowValues[j] = scaleValue*rowValues[j];
	  }
	  matrix_->replaceLocalValues( i, rowIndices, scaledRowValues );
	  rowIndices = Teuchos::null;
	  rowValues = Teuchos::null;
	}
      } else if( matrix_->getRangeMap()->isSameAs(*(x->getMap())) ) {
	return(-2); // Don't handle this case for now
      } else { // x.Map different than both matrix_->RowMap() and matrix_->ColMap()
	return(-2);
      }

      return(0);
    }

  //! @name Attribute access functions
  //@{

  //! Returns a pointer to the OperatorApplyInverse operator object that was used to create this InvOperator object.
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > getMatrix() const {return(matrix_);}

  //@}

 protected:
    Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > matrix_;
};

}

#endif /* TPETRA_SCALINGCRSMATRIX_H */
