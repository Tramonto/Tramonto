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

#ifndef DFT_2X2_A11_INVERTIBLE_SOLVER_H
#define DFT_2X2_A11_INVERTIBLE_SOLVER_H
class Epetra_Operator;
class Epetra_Vector;
class Epetra_Comm;
#include"dft_2x2_schur_epetra_operator.hpp"
#include "Epetra_Map.h"

//! dft_2x2_A11_invertible_solver:  Solver for a 2 by 2 block system with A11 explicitly invertible via offdiagonal term negation.
/*! Special solver for one important case of Tramonto.
  
*/
class dft_2x2_A11_invertible_solve{
    
  public:

  //@{ \name Constructors/destructor.
  //! dft_2x2_A11_invertible_solver Constructor.
  /* Initialize a solver for a 2 by 2 block system with A11 explicitly invertible via offdiagonal term negation.
     \param numBlocks (In) The number of physics blocks that will be handled by the solver manager.
     \param blockMatrix (In) 2-by-2 block of Epetra_CrsMatrix objects.
 */
  dft_2x2_A11_invertible_solver(int numBlocks, Epetra_Operator *** blockMatrix);

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
  updateMatrix(int numBlocks, Epetra_CrsMatrix *** blockMatrix)
  {op_.UpdateBlocks(blockMatrix[0][0], blockMatrix[0][1], blockMatrix[1][0], blockMatrix[1][1]);}
  //@}

  //@{ \name Apply solver method.
  //! Apply solver to given set of Rhs using initial guess in Lhs.
  /* Compute solution.
     \param blockRhs (In) 2 block of Epetra_Vector objects containing right hand side.
     \param blockLhs (InOut) 2 block of Epetra_Vector objects containing initial guess on entry, and solution on exit.
  */
  solve(Epetra_Vector ** blockRhs, Epetra_Vector ** blockLhs);
  //@}

private:
  dft_2x2_schur_epetra_operator op_;
  int numBlocks_;
  Epetra_CrsMatrix *** blockMatrix_;


#endif /* DFT_2X2_A11_INVERTIBLE_SOLVER_H */
