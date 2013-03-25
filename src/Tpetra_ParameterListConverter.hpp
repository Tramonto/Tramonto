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

#ifndef TPETRA_PARAMETERLISTCONVERTER_HPP
#define TPETRA_PARAMETERLISTCONVERTER_HPP

#include "BelosTypes.hpp"
#include "AztecOO.h"

namespace Tpetra {

  //! ParameterListConverter: A means of converting an Epetra ParameterList to a corresponding Tpetra ParameterList.

  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ParameterListConverter {
  public:

    /** \name Typedefs that give access to the template parameters. */
    //@{
    typedef typename Teuchos::ScalarTraits<Scalar>::halfPrecision halfScalar;
    typedef typename Teuchos::ScalarTraits<Scalar>::doublePrecision doubleScalar;

    //@}

    //! @name Constructor
    //@{
    //! Facilitates the conversion from Epetra parameters to Tpetra parameters.
    /*!
      \param In - A Teuchos ParameterList.
    */
    ParameterListConverter(Teuchos::RCP<Teuchos::ParameterList>& parameterListIn)
    {
      inputList_ = parameterListIn;
      isConverted_ = false;
      return;
    }

    //! Destructor
    virtual ~ParameterListConverter() {}
    //@}

    void convert()
    {
      Scalar sParam;
      //
      // Select scaling
      //
      int scaling  = inputList_->template get<int>( "Scaling" );
      outputList_.template set<int>( "Scaling",
				     scaling );
      //
      // Select preconditioner
      //
      int precond  = inputList_->template get<int>( "Precond" );
      double fill = inputList_->template get<double>("Ilut_fill");
      outputList_.template set<int>( "Precond",
				     precond );
      //
      // Ifpack2 Parameters
      //
      Teuchos::ParameterList ifpack2List;

      //
      // Use Ifpack2 Additive Schwarz
      //
      // Reordering
      int reorder = inputList_->template get<int>( "Reorder" );
      if (reorder == 1) {
	Teuchos::ParameterList zlist;
	zlist.set("order_method","rcm");
	ifpack2List.set( "schwarz: use reordering",
			 true );
	ifpack2List.set( "schwarz: reordering list",
			 zlist );
      }

      //
      // Use Ifpack2 ILUT on entire matrix or on subdomains with Additive Schwarz
      //
      // Level of fill
      ifpack2List.template set<Scalar>( "fact: ilut level-of-fill",
					inputList_->template get<double>("Ilut_fill") );
      // Absolute threshold
      ifpack2List.template set<Scalar>( "fact: absolute threshold",
					inputList_->template get<double>("Athresh") );
      // Relative threshold
      sParam = inputList_->template get<double>("Rthresh");
      if (sParam == 0.0)
	// The default value in AztecOO is 0.0; in Ifpack2 it is 1.0
	ifpack2List.template set<Scalar>( "fact: relative threshold", 1.0 );
      else
	ifpack2List.template set<Scalar>( "fact: relative threshold",
					  inputList_->template get<double>("Rthresh") );
      // Drop tolerance
      ifpack2List.template set<Scalar>( "fact: drop tolerance",
					inputList_->template get<double>("Drop") );

      outputList_.set( "ifpack2List", ifpack2List );

      //
      // Belos Parameters
      //
      Teuchos::ParameterList belosList;

      // Block size
      GlobalOrdinal blockSize = 1;
      belosList.set( "Block Size", blockSize);
      belosList.set( "Num Blocks",
		     inputList_->template get<int>("Kspace") / blockSize );
      // Maximum iterations
      belosList.set( "Maximum Iterations",
		     inputList_->template get<int>("Max_iter") );
      // Convergence tolerance
#if MIXED_PREC == 1
      belosList.template set<Scalar>( "Convergence Tolerance",
				      std::max( Teuchos::as<halfScalar>(inputList_->template get<double>("Tol")),
						5*Teuchos::ScalarTraits<halfScalar>::eps() ) );
#elif MIXED_PREC == 0
      belosList.template set<Scalar>( "Convergence Tolerance",
				      Teuchos::as<Scalar>(inputList_->template get<double>("Tol")) );
#endif

      belosList.set( "Orthogonalization",
		     "IMGS" );

      // Output
      //belosList.set( "Output Frequency", 10 );
      //belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
      //belosList.set( "Output Style", Belos::Brief);

      outputList_.set( "belosList", belosList );

      //
      // Muelu Parameters
      //
      Teuchos::ParameterList mueluList;

      // ML output
      mueluList.set( "ML output", 0 );

      // Smoother sweeps
      mueluList.set( "smoother: sweeps",  2 );

      // Smoother type
      mueluList.set( "smoother: type", "Chebyshev" );

      // Coarse sweeps
      mueluList.set( "coarse: sweeps",  6 );

      // Coarse type
      mueluList.set( "coarse: type",  "Chebyshev" );

      outputList_.set( "mueluList", mueluList );

      isConverted_ = true;
    }

    Teuchos::ParameterList getConvertedList()
    {
      TEUCHOS_TEST_FOR_EXCEPTION(isConverted_==false, std::runtime_error, "Parameter list has not be converted yet.\n");
      return outputList_;
    }

  protected:
    bool isConverted_;
    Teuchos::RCP<Teuchos::ParameterList> inputList_;
    Teuchos::ParameterList outputList_;
  };

}

#endif /* TPETRA_PARAMETERLISTCONVERTER_H */
