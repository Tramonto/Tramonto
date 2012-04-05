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

#ifndef DFT_POLYA22_TPETRA_BELOS_OPERATOR_H
#define DFT_POLYA22_TPETRA_BELOS_OPERATOR_H

#include "Tpetra_Headers.hpp"
#include "dft_PolyA22_Tpetra_Operator.hpp"
#include "dft_direct_solver_const.h"

//! dft_PolyA22_Tpetra_Operator: An implementation of the Operator class for Tramonto Schur complements with Coulomb effects.
/*! Special 2*numBeads by 2*numBeads (plus Coulomb) for Tramonto polymer problems.
*/    

template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal, 
         class Node = Kokkos::DefaultNode::DefaultNodeType>
class dft_PolyA22_Tpetra_Belos_Operator: 
  public virtual dft_PolyA22_Tpetra_Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> 
{
      
 public:
TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node);
  typedef dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> P22TO; 

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads (plus Coulomb) system

  dft_PolyA22_Tpetra_Belos_Operator
  (const RCP<const MAP> & cmsMap,const RCP<const MAP> & densityMap, 
   const RCP<const MAP> & block2Map, RCP<ParameterList> parameterList);

  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_PolyA22_Tpetra_Belos_Operator();
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a dft_PolyA22_Coulomb_Tpetra_Operator applied to a MultiVector X in Y.
    /*! 
    \param In
	   X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A MultiVector of dimension NumVectors containing result.
  */
  void 
  apply
  (const MV& X, MV& Y, Teuchos::ETransp mode = Teuchos::NO_TRANS, Scalar alpha = 1.0 , Scalar beta = 0.0) const;

  //! Returns the result of an inverse dft_PolyA22_Coulomb_Tpetra_Operator applied to a MultiVector X in Y.
  /*! 
    \param In
    X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y - A MultiVector of dimension NumVectors containing result.
  */
  void
  applyInverse
  (const MV& X,MV& Y) const;

  //@}

};
#endif /* DFT_POLYA22_TPETRA_BELOS_OPERATOR_H */
