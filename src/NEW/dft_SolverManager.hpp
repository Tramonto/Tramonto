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
#ifndef DFT_SOLVERMANAGER_HPP
#define DFT_SOLVERMANAGER_HPP

class Epetra_Operator;
class Epetra_Vector;
class Epetra_Comm;
#include "Epetra_Map.h"

//! dft_SolverManager:  Solver manager class for Tramonto using Trilinos.
/*! The dft_SolverManager class supports solver capabilities for Tramonto.
  
*/
class dft_SolverManager {
    
  public:
  //@{ \name Constructors/destructor.
  //! dft_SolverManager Constructor.
  /* Initialize a solver manager for Tramonto
     \param numBlocks (In) The number of physics blocks that will be handled by the solver manager.
     \param comm (In) Epetra communicator that should be used by the solver.
  */
  dft_SolverManager(int numBlocks, const Epetra_Comm & comm);

  //! dft_SolverManager Destructor.
  /*! Completely deletes a dft_SolverManager object.
  */
  virtual ~dft_SolverManager();
  //@}

  //@{ \name Block structure setup methods

  //! Define global to local index row mappings
  /* Define the list of global IDs for the ith physics type.  This mapping is used for identifying
     the global row indexing given the local row IDs.
     \param i (In) Physics type block number
     \param numGIDs (In) Number of global IDs
     \param GIDs (In) List of global IDs in the order of local indexing.
  */
  int setRowMap(int i, int numGIDs int * GIDs) {rowMap_[i] = new Epetra_Map(-1, numGIDs, GIDs, 0, comm_); return;};

  //! Define global to local index column mappings
  /* Define the list of global IDs for the jth physics type.  This mapping is used for identifying
     the global column indexing given the local column IDs.
     \param i (In) Physics type block number
     \param numGIDs (In) Number of global IDs
     \param GIDs (In) List of global IDs in the order of local indexing.
  */
  int setColMap(int j, int numGIDs int * GIDs) {colMap_[j] = new Epetra_Map(-1, numGIDs, GIDs, 0, comm_); return;};
  //! Method that must be called once, when all row and column maps are set.
  int finalizeBlockStructure();
  //@}

  //@{ \name Graph setup methods

  //! Insert graph indices
  /* Insert indices into block i,j.
     \param i (In) ith physics block.
     \param j (In) jth physics block (column).
     \param localRow (In) Local row to which entries will be added.  
            This value should be between zero and the number of equations in block i.
     \param numEntries (In) Number of graph indices being inserted.
     \param indices (In) Column indices corresponding to matrix values.
  */
  int insertRowIndices(int i, int j, int localRow, int numEntries, int * indices) {
    EPETRA_CHK_ERR(blockGraph_[i][j]->InsertMyIndices(localRow, numEntries, indices));
    return(0);
  }

  //! Method that must be called once, when all graph indices are set.
  int finalizeGraphStructure();
  //@}

  //@{ \name Matrix value setup methods

  //! Method that must be called each time \e prior to starting matrix value insertion (usually once per nonlinear iteration).
  int initializeMatrixFill();

  //! Insert matrix coefficients
  /* Insert values into block i,j.
     \param i (In) ith physics block.
     \param j (In) jth physics block (column).
     \param localRow (In) Local row to which entries will be added.  
            This value should be between zero and the number of equations in block i.
     \param numEntries (In) Number of matrix coefficients being inserted.
     \param values (In) Matrix values.
     \param indices (In) Column indices corresponding to matrix values.
  */
  int insertRowValues(int i, int j, int localRow, int numEntries, double * values, int * indices) {
    EPETRA_CHK_ERR(blockMatrix_[i][j]->SumIntoMyValues(localRow, numEntries, values, indices));
    return(0);
  }
  //! Method that must be called each time matrix value insertion is complete (usually once per nonlinear iteration).
  int finalizeMatrixFill();
  //@}

  //@{ \name Solver Strategy Methods
  //! Select a strategy for solving the problem;  for now we just hard-code the strategy.
  /*! This method is typically called once, or a few times at most per run of the application.
      All of the interesting future work will go into developing good strategies using the solver
      tools we have on hard, or are actively developing.
  */
  int setSolverStrategy(){return(0)} // This is a no-op for now.

  //! Setup up the solver for solving the current linear problem.
  int setupSolver();
  //@}

private:

  Epetra_Comm comm_;
  Epetra_Matrix *** blockMatrix_;
  Epetra_Graph *** blockGraph_;
  Epetra_Vector ** blockLhs_;
  Epetra_Vector ** blockRhs_;
  Epetra_Map ** rowMaps_;
  Epetra_Map ** colMaps_;
  bool isStructureSet_;

};

#endif /* DFT_SOLVERMANAGER_HPP */
