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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ***********************************************************************
//@HEADER
*/
#ifndef DFT_2X2_A11_INVERTIBLE_SOLVER_H
#define DFT_2X2_A11_INVERTIBLE_SOLVER_H
class Epetra_Operator;
class Epetra_Vector;
class Epetra_Comm;
#include "Epetra_Map.h"

//! dft_2x2_A11_invertible_solver:  Solver for a 2 by 2 block system with A11 explicitly invertible via offdiagonal term negation.
/*! Special solver for one important case of Tramonto.
  
*/
class dft_2x2_A11_invertible_solver {
    
  public:

  //@{ \name Constructors/destructor.
  //! dft_2x2_A11_invertible_solver Constructor.
  /* Initialize a solver for a 2 by 2 block system with A11 explicitly invertible via offdiagonal term negation.
     \param numBlocks (In) The number of physics blocks that will be handled by the solver manager.
     \param blockMatrix (In) 2-by-2 block of Epetra_CrsMatrix objects.
 */
  dft_2x2_A11_invertible_solver(int numBlocks, Epetra_CrsMatrix *** blockMatrix);

  //! dft_2x2_A11_invertible_solver Destructor.
  /*! Completely deletes a dft_2x2_A11_invertible_solver object.
  */
  virtual ~dft_2x2_A11_invertible_solver();
  //@}

  //@{ \name Update method.
  //! Update dft_2x2_A11_invertible_solver with new values (same structure).
  /* Update a solver for a 2 by 2 block system with A11 explicitly invertible via offdiagonal term negation.
     \param numBlocks (In) The number of physics blocks that will be handled by the solver manager.
     \param blockMatrix (In) 2-by-2 block of Epetra_CrsMatrix objects.
  */
  updateMatrix(int numBlocks, Epetra_CrsMatrix *** blockMatrix);
  //@}

  //@{ \name Apply solver method.
  //! Apply solver to given set of Rhs using initial guess in Lhs.
  /* Compute solution.
     \param blockRhs (In) 2 block of Epetra_Vector objects containing right hand side.
     \param blockLhs (InOut) 2 block of Epetra_Vector objects containing initial guess on entry, and solution on exit.
  */
  apply(Epetra_Vector ** blockLhs, Epetra_Vector ** blockRhs);
  //@}


#endif /* DFT_2X2_A11_INVERTIBLE_SOLVER_H */
