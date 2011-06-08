//@HEADER
// ********************************************************************
// Tramonto: A molecular theory code for structured and uniform fluids
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation; either version 2.1
// of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// 02110-1301, USA.
// ********************************************************************
//@HEADER

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
  OperatorApplyInverse<Scalar,LocalOrdinal,GlobalOrdinal,Node>::~OperatorApplyInverse() {
  }


} // namespace Tpetra


#endif //TPETRA_OPERATORAPPLYINVERSE_HPP
