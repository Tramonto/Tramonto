/*@HEADER
// ***********************************************************************
// 
//                Tramonto: Molecular Theories Modeling Code
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// Questions? Contact Laura J.D. Frink (ljfrink@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/
#ifndef DFT_LINSOLVE_H
#define DFT_LINSOLVE_H

#include "dft_vec_space.h"

/*! \file dft_linsolve.h
\brief Function prototypes required for Tramonto linear solver functionality

The file must be included by each Tramonto source file that needs access
to linear solver functionality.  This file also specifies the functions that a 
linear solver library must implement in order to supply solver functionality
for Tramonto.

*/
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

/*! \fn int dft_linsolve_set_prec_options(DFT_LINSOLVE_MGR * mgr, int option_key, int option_value)

\brief Set an option value for the solver manager to use with a preconditioner.

\param mgr (In) A pointer to a valid manager object.
\param option_key (In) Option key.
\param option_value (In) Value to associate with specified key.
\return Return a 0 if option key/value pair was understood and accepted by manager, -1 if not understood, +1 if understood but not accepted.
*/

void dft_linsolve_set_prec_options(DFT_LINSOLVE_MGR * mgr, int option_key, int option_value);


#endif /* DFT_LINSOLVE_H */
