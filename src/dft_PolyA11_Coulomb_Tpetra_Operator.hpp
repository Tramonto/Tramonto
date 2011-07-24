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

#ifndef DFT_POLYA11_COULOMB_TPETRA_OPERATOR_H
#define DFT_POLYA11_COULOMB_TPETRA_OPERATOR_H

#include "Tpetra_Headers.hpp"
#include "dft_PolyA11_Tpetra_Operator.hpp"
#include "dft_direct_solver_const.h"

//! dft_PolyA11_Coulomb_Tpetra_Operator: An implementation of the Operator class for Tramonto Schur complements with Coulomb effects.
/*! Special 2*numBeads by 2*numBeads (plus Coulomb) for Tramonto polymer problems.
*/    

template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal, 
         class Node = Kokkos::DefaultNode::DefaultNodeType>
class dft_PolyA11_Coulomb_Tpetra_Operator: 
  public virtual dft_PolyA11_Tpetra_Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> 
{
      
 public:
TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node);
  typedef dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> P11TO;

  //@{ \name Constructors.
    //! Builds an implicit composite operator from a 2*numBeads by 2*numBeads (plus Coulomb) system
  /*  dft_PolyA11_Coulomb_Operator(const Map & ownedMap, const Map & block1Map, const Map & allGMap, const Map & poissonMap, LocalOrdinal * solverOptions, Scalar * solverParams);*/

  dft_PolyA11_Coulomb_Tpetra_Operator
  (RCP<const MAP > & ownedMap,  RCP<const MAP > & block1Map, RCP<const MAP > & allGMap, 
   RCP<const MAP > & poissonMap, RCP<ParameterList> parameterList);
  //@}
  //@{ \name Assembly methods.
  void
  initializeProblemValues
  ();
 
  void
  insertMatrixValue
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, GlobalOrdinal rowGID, GlobalOrdinal colGID, Scalar value);

  void
  finalizeProblemValues
  ();
  //@}
  //@{ \name Destructor.
    //! Destructor
  virtual ~dft_PolyA11_Coulomb_Tpetra_Operator();
  //@}

  //@{ \name Mathematical functions.

    //! Returns the result of a dft_PolyA11_Coulomb_Tpetra_Operator applied to a MultiVector X in Y.
    /*! 
    \param In
	   X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A MultiVector of dimension NumVectors containing result.
  */
  void 
  apply
  (const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const;

  //! Returns the result of an inverse dft_PolyA11_Coulomb_Tpetra_Operator applied to a MultiVector X in Y.
  /*! 
    \param In
    X - A MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
    Y - A MultiVector of dimension NumVectors containing result.
  */
  void
  applyInverse
  (const MV& X, MV& Y) const;

  //@}
  
  //@{ \name Atribute access functions
  //! Returns a pointer to the Comm communicator associated with this operator.
  const RCP<const COMM> & 
  Comm
  () const
  {
    return(block1Map_->getComm());
  };
  
  //! Returns the Map object associated with the domain of this operator.
  const RCP<const MAP> & 
  getDomainMap
  () const 
  { 
    return(block1Map_);
  };
  
  //! Returns the Map object associated with the range of this operator.
  const RCP<const MAP> & 
  getRangeMap
  () const 
  {
    return(block1Map_);
  };
  //@}

protected:

  void insertPoissonRow();
  const RCP<const MAP> block1Map_;
  const RCP<const MAP> allGMap_;
  const RCP<const MAP> poissonMap_;
  RCP<ParameterList> parameterList_;
  RCP<MAT> poissonMatrix_;
  GlobalOrdinal curPoissonRow_;
  GlobalOrdinal curPoissonOwnedNode_;
  std::map<GlobalOrdinal, Scalar> curPoissonRowValues_;
  using P11TO::isGraphStructureSet_;
  using P11TO::Label_;
  using P11TO::isLinearProblemSet_;
  using P11TO::firstTime_;
  using P11TO::numBlocks_;
  using P11TO::curRow_;
  using P11TO::matrix_;
  using P11TO::curOwnedPhysicsID_;
  using P11TO::curRowValues_;
  using P11TO::indices_;
  using P11TO::values_;
  using P11TO::ownedMap_;
  using P11TO::curOwnedNode_;

};
#endif /* DFT_POLYA11_COULOMB_TPETRA_OPERATOR_H */
