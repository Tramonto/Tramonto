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
#ifndef DFT_POLYLINPROBMGR_HPP
#define DFT_POLYLINPROBMGR_HPP

class Epetra_Vector;
class Epetra_Map;
class Epetra_Import;
class Epetra_CrsMatrix;
class Epetra_LinearProblem;
class AztecOO;
#include "dft_BasicLinProbMgr.hpp"
#include "Epetra_MpiComm.h"
#include "Teuchos_RefCountPtr.hpp"

//! dft_PolyLinProbMgr:  Solver manager class for Tramonto using Trilinos.
/*! The dft_PolyLinProbMgr class supports solver capabilities for Tramonto.
  
*/
class dft_PolyLinProbMgr: public dft_BasicLinProbMgr {
    
  public:
  //@{ \name Constructors/destructor.
  //! dft_PolyLinProbMgr Constructor.
  /*! Initialize a linear problem manager for Tramonto
     \param numUnknownsPerNode (In) The number of unknowns tracked per node of the mesh.
     \param solverOptions (In) An array of ints defined in dft_solver_defs.h containing information to 
                               guide and report solver status.
     \param solverParams (In) An array of doubles defined in dft_solver_defs.h containing information to 
                               guide and report solver status.
     \param comm (In) MPI communicator that should be used by the solver.
  */
  dft_PolyLinProbMgr(int numUnknownsPerNode, int * solverOptions, double * solverParams, MPI_Comm comm);

  //! dft_PolyLinProbMgr Destructor.
  /*! Completely deletes a dft_PolyLinProbMgr object.
  */
  virtual ~dft_PolyLinProbMgr();
  //@}

  //@{ \name Block structure setup methods

  //! Define G equation IDs
  /*! Define the list of physics IDs associated with the G equations in ascending order..
     \param numGEquations (In) Number of G equations per node.
     \param physicsIDs (In) List of physics IDs associated with the G equations in ascending order.
  */
  int setGEquationIDs(int numGEquations, int * physicsIDs) { 
    gEquations_.Resize(numGEquations); 
    for (int i=0; i<numGEquations; i++) gEquations_[i] = physicsIDs[i];
  }

  //! Define G inverse equation IDs
  /*! Define the list of physics IDs associated with the G inverse equations in descending order.
     \param numGInvEquations (In) Number of G inverse equations per node.
     \param physicsIDs (In) List of physics IDs associated with the G equations in ascending order.
  */
  int setGInvEquationIDs(int numGInvEquations, int * physicsIDs) { 
    gInvEquations_.Resize(numGInvEquations); 
    for (int i=0; i<numGInvEquations; i++) gInvEquations_[i] = physicsIDs[i];
  }

  //! Define CMS equation IDs
  /*! Define the list of physics IDs associated with the CMS equations in ascending order.
     \param numCmsEquations (In) Number of CMS equations per node.
     \param physicsIDs (In) List of physics IDs associated with the CMS equations in ascending order.
  */
  int setCmsEquationIDs(int numCmsEquations, int * physicsIDs) { 
    cmsEquations_.Resize(numCmsEquations); 
    for (int i=0; i<numCmsEquations; i++) cmsEquations_[i] = physicsIDs[i];
  }

  //! Define primitive density equation IDs
  /*! Define the list of physics IDs associated with the primitive density equations in ascending order.
     \param numDensityEquations (In) Number of primitive density equations per node.
     \param physicsIDs (In) List of physics IDs associated with the primitive density equations in ascending order.
  */
  int setDensityEquationIDs(int numDensityEquations, int * physicsIDs) { 
    densityEquations_.Resize(numDensityEquations); 
    for (int i=0; i<numDensityEquations; i++) densityEquations_[i] = physicsIDs[i];
  }

  //! Method that must be called once, when all row and column maps are set.
  /*! This method constructs all of the Epetra_CrsGraph objects and the lhs and rhs vectors. */
  int finalizeBlockStructure();
  //@}

  //@{ \name Matrix, lhs and rhs value setup methods for piece-wise construction.

  //! Method that must be called each time \e prior to starting matrix, lhs and rhs value insertion (usually once per nonlinear iteration).
  /*! This method zeros out the matrix, lhs and rhs values. */
  virtual int initializeProblemValues();

  //! Insert rhs value based on ownedNode and ownedPhysicsID.
  /*! Insert rhs value into entry based on owned node and physicsID.
     \param ownedPhysicsID (In) The index for the type of unknown.  
                           This should be between 0 and one less than the number of physics variables tracked at a node.  
     \param ownedNode (In) Current owned node ID.  This is the local ID based on the set of owned nodes for this processor.
                      This should be between 0 and one less than the number of nodes owned by this processor, independent of
		      the number of physics types being computed for the given problem.
     \param value (In) Rhs value.
  */
  virtual int insertRhsValue(int ownedPhysicsID, int ownedNode, double value);

  //! Insert single matrix coefficient into system.
  /*! Insert single value into matrix.
     \param ownedPhysicsID (In) The physics ID for the matrix row being updated.  
                           This should be between 0 and one less than the number of physics variables tracked at a node.  
			   ownedNode and ownedPhysicsID together specify which row of the matrix is being updated.
     \param ownedNode (In) Current owned node ID.  This is the local ID based on the set of owned nodes for this processor.
                      This should be between 0 and one less than the number of nodes owned by this processor, independent of
		      the number of physics types being computed for the given problem. rowPhysicsID and ownedNode together specify 
		      which row of the matrix is being updated.
     \param boxPhysicsID (In) The index for the type of unknown to use for the column indices.  
                           This should be between 0 and one less than the number of physics variables tracked at a node.  
			   boxNode and colPhysicsID together specify which column of the matrix is being updated.
     \param boxNode (In) Current box node ID.  Local ID based on the box node ordering.
     \param value (In) Matrix value.
  */
  virtual int insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value);

  //! Method that must be called each time matrix value insertion is complete (usually once per nonlinear iteration).
  virtual int finalizeProblemValues();

  //@}


  //@{ \name Solver-related methods

  //! Setup up the solver for solving the current linear problem.
  virtual int setupSolver();

  //! Solve the current linear problem.
  virtual int solve();
  //@}

  //@{ \name Output facilities

  //! Write matrix to specified filename using Matrix Market (i,j,value) format.
  /*! Write the matrix owned by the solver manager to the file called filename on PE 0 regardless of how many processors are
      involved in the computation.
      \param filename (In) Name of the file.  Any contents in this file will be erased.
      \param matrixName (In) Optional brief name for the matrix, will be inserted in the header of the output file, can be set to 0.
      \param matrixDescription (In) Optional longer description for the matrix, will be inserted in the header of the output file, can be set to 0.

    \return Returns 0 if no error, -1 of any problems with the file system.
  */
  virtual int writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const;
  //@}

  //@{ \name Miscellaneous support methods (used by the application, or by LOCA, or both)

  //! Apply global linear operator for all physics types, b = Ax
  /*! Computes b = Ax where A is the global linear operator and each x[i] and b[i] is associated with the ith physics block.
    \param x (In) Array of pointers such that x[i] is an array of doubles of length numBoxNodes.
    \param b (Out) Array of pointers such that b[i] is an array of doubles of length numOwnedNodes.

    \warning Note that x[i] is of length numBoxNodes and b[i] is of length numOwnedNodes.
  */
  virtual int applyMatrix(const double** x, double** b) const;
  
  //@}

protected:

  Epetra_IntSerialDenseVector gEquations_;
  Epetra_IntSerialDenseVector gInvEquations_;
  Epetra_IntSerialDenseVector cmsEquations_;
  Epetra_IntSerialDenseVector densityEquations_;
  Epetra_IntSerialDenseVector physicsOrdering_;
  Epetra_IntSerialDenseVector physicsIdToMatrixId_;
  Epetra_IntSerialDenseVector physicsIdToSchurBlockId_;
  
	     
};

#endif /* DFT_POLYLINPROBMGR_HPP */
