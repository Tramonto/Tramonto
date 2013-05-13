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
#include "az_aztec_defs.h"

namespace Tpetra {

  //! ParameterListConverter: A means of converting an Epetra ParameterList to a corresponding Tpetra ParameterList.

  template<class Scalar, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
  class ParameterListConverter {
  public:

    /** \name Typedefs that give access to the template parameters. */
    //@{
    typedef typename Teuchos::ScalarTraits<Scalar>::halfPrecision halfScalar;
    typedef typename Teuchos::ScalarTraits<Scalar>::doublePrecision doubleScalar;
    typedef typename Teuchos::ScalarTraits<Scalar> STS;
    typedef typename Teuchos::ScalarTraits<halfScalar> STHS;

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
      outputList_.template set<int>( "Precond",
				     inputList_->template get<int>( "Precond" ) );
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
      ifpack2List.template set<double>( "fact: ilut level-of-fill",
      					inputList_->template get<double>("Ilut_fill") );

      // Absolute threshold
      ifpack2List.template set<double>( "fact: absolute threshold",
					inputList_->template get<double>("Athresh") );
      // Relative threshold
      sParam = inputList_->template get<double>("Rthresh");
      if (sParam == 0.0)
	// The default value in AztecOO is 0.0; in Ifpack2 it is 1.0
	ifpack2List.template set<double>( "fact: relative threshold", 1.0 );
      else
	ifpack2List.template set<double>( "fact: relative threshold",
					  inputList_->template get<double>("Rthresh") );
      // Drop tolerance
      sParam = inputList_->template get<double>("Drop");
      if (sParam == 0.0) {
	// Don't set the drop tolerence, let Ifpack2 choose the default
      } else {
	ifpack2List.template set<double>( "fact: drop tolerance",
					  inputList_->template get<double>("Drop") );
      }

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
				      Teuchos::as<halfScalar>(inputList_->template get<double>("Tol")));
#else
      belosList.template set<Scalar>( "Convergence Tolerance",
				      Teuchos::as<Scalar>(inputList_->template get<double>("Tol")) );
#endif

      // Orthogonalization
      int orthog  = inputList_->template get<int>( "Orthog" );
      if (orthog == AZ_modified) {
	// Modified Gram Schmidt
	belosList.set( "Orthogonalization",
		       "IMGS" );
      } else {
	// Classic Gram Schmidt
	belosList.set( "Orthogonalization",
		       "ICGS" );
      }

      // Output
      // belosList.set( "Output Frequency", 20 );
      // belosList.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
      // belosList.set( "Output Style", Belos::Brief);

      outputList_.set( "belosList", belosList );

      Teuchos::ParameterList belosListSchur;

      // Block size
      GlobalOrdinal blockSizeSchur = 1;
      GlobalOrdinal maxIterSchur = 500;
      belosListSchur.set( "Block Size", blockSizeSchur);
      belosListSchur.set( "Num Blocks",
			  maxIterSchur / blockSizeSchur );
      // Maximum iterations
      belosListSchur.set( "Maximum Iterations",
			  maxIterSchur );

      // Convergence tolerance
#if MIXED_PREC == 1
      belosListSchur.template set<Scalar>( "Convergence Tolerance",
					   Teuchos::as<halfScalar>(inputList_->template get<double>("Tol")));
#else
      belosListSchur.template set<Scalar>( "Convergence Tolerance",
					   Teuchos::as<Scalar>(inputList_->template get<double>("Tol")) );
#endif

      // Orthogonalization
      int orthogSchur  = inputList_->template get<int>( "Orthog" );
      if (orthogSchur == AZ_modified) {
	// Modified Gram Schmidt
	belosListSchur.set( "Orthogonalization",
			    "IMGS" );
      } else {
	// Classic Gram Schmidt
	belosListSchur.set( "Orthogonalization",
			    "ICGS" );
      }

      // Output
      //belosListSchur.set( "Output Frequency", 10 );
      //belosListSchur.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
      //belosListSchur.set( "Output Style", Belos::Brief);

      outputList_.set( "belosListSchur", belosListSchur );

      Teuchos::ParameterList belosListA11;

      // Block size
      GlobalOrdinal blockSizeA11 = 1;
      GlobalOrdinal maxIterA11 = 500;
      belosListA11.set( "Block Size", blockSizeA11);
      belosListA11.set( "Num Blocks",
			maxIterA11 / blockSizeA11 );
      // Maximum iterations
      belosListA11.set( "Maximum Iterations",
			maxIterA11 );

      // Convergence tolerance
#if MIXED_PREC == 1
      belosListA11.template set<Scalar>( "Convergence Tolerance",
					 Teuchos::as<halfScalar>(inputList_->template get<double>("Tol")));
#else
      belosListA11.template set<Scalar>( "Convergence Tolerance",
					 Teuchos::as<Scalar>(inputList_->template get<double>("Tol")) );
#endif

      // Orthogonalization
      int orthogA11  = inputList_->template get<int>( "Orthog" );
      if (orthogA11 == AZ_modified) {
	// Modified Gram Schmidt
	belosListA11.set( "Orthogonalization",
			  "IMGS" );
      } else {
	// Classic Gram Schmidt
	belosListA11.set( "Orthogonalization",
			  "ICGS" );
      }

      // Output
      //belosListA11.set( "Output Frequency", 10 );
      //belosListA11.set( "Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails );
      //belosListA11.set( "Output Style", Belos::Brief);

      outputList_.set( "belosListA11", belosListA11 );

      //
      // Muelu Parameters
      //
      Teuchos::ParameterList fillCompleteList;

      // Promise that no nonlocal changes have been made
      // This prevents unnecessary communication in fillComplete()
      // fillCompleteList.set( "No Nonlocal Changes", true );

      outputList_.set( "fillCompleteList", fillCompleteList );

      //
      // Muelu Parameters
      //
      Teuchos::ParameterList mueluList;

      // ML output
      //      mueluList.set( "ML output", 0 );

      // Smoother sweeps
      mueluList.set( "smoother: sweeps",  2 );

      // Smoother type
      mueluList.set( "smoother: type", "Chebyshev" );

      // Coarse sweeps
      mueluList.set( "coarse: sweeps",  6 );

      // Coarse type
      mueluList.set( "coarse: type",  "Chebyshev" );

      outputList_.set( "mueluList", mueluList );

      //
      // Use Ifpack2 ILUT on subdomains with Additive Schwarz for A22 block
      //
      Teuchos::ParameterList ifpack2ListA22;
      // Reordering
      int reorderA22 = inputList_->template get<int>( "Reorder" );
      if (reorderA22 == 1) {
	Teuchos::ParameterList zlist;
	zlist.set("order_method","rcm");
	ifpack2ListA22.set( "schwarz: use reordering",
			    true );
	ifpack2ListA22.set( "schwarz: reordering list",
			    zlist );
      }

      // Level of fill
      ifpack2ListA22.template set<double>( "fact: ilut level-of-fill",
					   7.0 );

      // Absolute threshold
      ifpack2ListA22.template set<double>( "fact: absolute threshold",
					   inputList_->template get<double>("Athresh") );
      // Relative threshold
      sParam = inputList_->template get<double>("Rthresh");
      if (sParam == 0.0)
	// The default value in AztecOO is 0.0; in Ifpack2 it is 1.0
	ifpack2ListA22.template set<double>( "fact: relative threshold", 1.0 );
      else
	ifpack2ListA22.template set<double>( "fact: relative threshold",
					     inputList_->template get<double>("Rthresh") );
      // Drop tolerance
      sParam = inputList_->template get<double>("Drop");
      if (sParam == 0.0) {
	// Don't set the drop tolerence, let Ifpack2 choose the default
      } else {
	ifpack2ListA22.template set<double>( "fact: drop tolerance",
					     inputList_->template get<double>("Drop") );
      }

      outputList_.set( "ifpack2ListA22", ifpack2ListA22 );

      //
      // Use Ifpack2 ILUT on subdomains with Additive Schwarz for A11 block
      //
      Teuchos::ParameterList ifpack2ListA11;

      // Level of fill
      ifpack2ListA11.template set<double>( "fact: ilut level-of-fill",
					   inputList_->template get<double>("Ilut_fill") );

      // Absolute threshold
      ifpack2ListA11.template set<double>( "fact: absolute threshold",
					   inputList_->template get<double>("Athresh") );
      // Relative threshold
      sParam = inputList_->template get<double>("Rthresh");
      if (sParam == 0.0)
	// The default value in AztecOO is 0.0; in Ifpack2 it is 1.0
	ifpack2ListA11.template set<double>( "fact: relative threshold", 1.0 );
      else
	ifpack2ListA11.template set<double>( "fact: relative threshold",
					     inputList_->template get<double>("Rthresh") );
      // Drop tolerance
      sParam = inputList_->template get<double>("Drop");
      if (sParam == 0.0) {
	// Don't set the drop tolerence, let Ifpack2 choose the default
      } else {
	ifpack2ListA11.template set<double>( "fact: drop tolerance",
					     inputList_->template get<double>("Drop") );
      }

      outputList_.set( "ifpack2ListA11", ifpack2ListA11 );

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
