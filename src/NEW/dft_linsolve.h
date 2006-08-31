/*
//@HEADER
// ********************************************************************
// Copyright (2006) Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000, there is a non-exclusive license for use of this
// work by or on behalf of the U.S. Government. Export of this program
// may require a license from the United States Government.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// ********************************************************************
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
