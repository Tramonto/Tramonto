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

#include "dft_PolyA22_Coulomb_Tpetra_Belos_Operator.hpp"

//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Coulomb_Tpetra_Belos_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_PolyA22_Coulomb_Tpetra_Belos_Operator
(const RCP<const MAP> & cmsMap, const RCP<const MAP> & densityMap,
 const RCP<const MAP> & poissonMap, const RCP<const MAP> & cmsDensMap,
 const RCP<const MAP> & block2Map,  RCP<ParameterList> parameterList) :
  dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  (cmsMap, densityMap, cmsDensMap, parameterList),
  dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>
  (cmsMap, densityMap, poissonMap, cmsDensMap, block2Map, parameterList)
{
} //end constructor
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_PolyA22_Coulomb_Tpetra_Belos_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_PolyA22_Coulomb_Tpetra_Belos_Operator
()
{
  return;
} //end destructor
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Belos_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyInverse
(const MV& X, MV& Y) const
{
  P22CO::apply( X, Y );
} //end applyInverse
//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_PolyA22_Coulomb_Tpetra_Belos_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
apply
(const MV& X, MV& Y, Teuchos::ETransp mode, Scalar alpha, Scalar beta) const
{
  P22CO::applyInverse( X, Y );
} //end Apply
//==============================================================================

#if LINSOLVE_PREC == 0
// Use float
template class dft_PolyA22_Coulomb_Tpetra_Belos_Operator<float, int, int>;
#elif LINSOLVE_PREC == 1
// Use double
template class dft_PolyA22_Coulomb_Tpetra_Belos_Operator<double, int, int>;
#elif LINSOLVE_PREC == 2
// Use quad double
template class dft_PolyA22_Coulomb_Tpetra_Belos_Operator<qd_real, int, int>;
#elif LINSOLVE_PREC == 3
// Use double double
template class dft_PolyA22_Coulomb_Tpetra_Belos_Operator<dd_real, int, int>;
#endif
