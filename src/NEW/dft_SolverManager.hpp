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
#ifndef DFT_SOLVERMANAGER_HPP
#define DFT_SOLVERMANAGER_HPP

//#include "dft_solver_defs.h"
class Epetra_MpiComm;
class Epetra_RowMatrix;
class Epetra_CrsGraph;
class Epetra_Vector;
class Epetra_BlockMap;
class Epetra_Map;
class Epetra_Import;
class Epetra_CrsMatrix;
class AztecOO;
#include "Teuchos_RefCountPtr.hpp"

//! dft_SolverManager:  Solver manager class for Tramonto using Trilinos.
/*! The dft_SolverManager class supports solver capabilities for Tramonto.
  
*/
class dft_SolverManager {
    
  public:
  //@{ \name Constructors/destructor.
  //! dft_SolverManager Constructor.
  /*! Initialize a solver manager for Tramonto
     \param numUnknownsPerNode (In) The number of unknowns tracked per node of the mesh.
     \param unknownToPhysicsType (In) An array of length numUnknownsPerNode, such that unknownToPhysicsType[i] 
            describes the ith physics type.
     \param solverOptions (In) An array of ints defined in dft_solver_defs.h containing information to 
                               guide and report solver status.
     \param solverParams (In) An array of doubles defined in dft_solver_defs.h containing information to 
                               guide and report solver status.
     \param comm (In) MPI communicator that should be used by the solver.
  */
  dft_SolverManager(int numUnknownsPerNode, int * unknownToPhysicsType, int * solverOptions, double * solverParams, MPI_Comm comm);

  //! dft_SolverManager Destructor.
  /*! Completely deletes a dft_SolverManager object.
  */
  virtual ~dft_SolverManager();
  //@}

  //@{ \name Block structure setup methods

  //! Define global to local index row mappings
  /*! Define the list of global node IDs owned by this processor.  This mapping is used for identifying
     the global row indexing given the local row IDs and the physics ID.
     \param numOwnedNodes (In) Number of owned nodes.
     \param GIDs (In) List of global IDs in the order of local indexing of the owned nodes.
     \param nx (In) Number of nodes in X direction on this processor.  If set to 0 (default value), nx, ny and nz 
                    will be ignored.  For a 1D problem, this value should be set to numGIDs.
     \param ny (In) Number of nodes in Y direction on this processor, defaults to 1.
     \param nz (In) Number of nodes in Z direction on this processor, defaults to 1.

  */
  int setNodalRowMap(int numOwnedNodes, int * GIDs, int nx=0, int ny = 1, int nz = 1);

  //! Define global to local index column mappings, the rectangular box containing all ghost nodes and owned nodes.
  /*! Define the list of global node IDs reached to by global row nodes on this processor. 
     This mapping is used for identifying
     the global column indexing given the local column IDs and the physics ID.
     \param numBoxNodes (In) Number of global IDs.  Should equal nx*ny*nz.
     \param GIDs (In) List of global IDs in the order of local indexing.
     \param nx (In) Number of nodes in X direction on this processor.  
                    For a 1D problem or unstructured problem, this value should be set to numGIDs.
     \param ny (In) Number of nodes in Y direction on this processor, defaults to 1.
     \param nz (In) Number of nodes in Z direction on this processor, defaults to 1.

  */
  int setNodalColMap(int numBoxNodes, int * GIDs, int nx=0, int ny = 1, int nz = 1);

  //! Method that must be called once, when all row and column maps are set.
  /*! This method constructs all of the Epetra_CrsGraph objects and the lhs and rhs vectors. */
  int finalizeBlockStructure();
  //@}

  //@{ \name Matrix, lhs and rhs value setup methods for piece-wise construction.

  //! Method that must be called each time \e prior to starting matrix, lhs and rhs value insertion (usually once per nonlinear iteration).
  /*! This method zeros out the matrix, lhs and rhs values. */
  int initializeProblemValues();

  //! Insert rhs value based on ownedNode and ownedPhysicsID.
  /*! Insert rhs value into entry based on owned node and physicsID.
     \param ownedPhysicsID (In) The index for the type of unknown.  
                           This should be between 0 and one less than the number of physics variables tracked at a node.  
     \param ownedNode (In) Current owned node ID.  This is the local ID based on the set of owned nodes for this processor.
                      This should be between 0 and one less than the number of nodes owned by this processor, independent of
		      the number of physics types being computed for the given problem.
     \param value (In) Rhs value.
  */
  int insertRhsValue(int ownedPhysicsID, int ownedNode, double value);

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
  int insertMatrixValue(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int boxNode, double value);

  //! Insert matrix coefficients for a given row, where columns are all from the same physics type at different nodes.
  /*! Insert values into matrix for a given row, where columns are all from the same physics type at different nodes.
     \param ownedPhysicsID (In) The physics ID for the matrix row being updated.
                           This should be between 0 and one less than the number of physics variables tracked at a node.  
			   ownedNode and rowPhysicsID together specify which row of the matrix is being updated.
     \param ownedNode (In) Current owned node ID.  This is the local ID based on the set of owned nodes for this processor.
                      This should be between 0 and one less than the number of nodes owned by this processor, independent of
		      the number of physics types being computed for the given problem.
     \param boxPhysicsID (In) The index for the type of unknown to use for the column indices.  
                           This should be between 0 and one less than the number of physics variables tracked at a node.  
			   boxNodelist and colPhysicsID together specify which column of the matrix is being updated.
     \param boxNodeList (In) Nodal Column indices corresponding to matrix values. These are local IDs in the box coordinate system.
     \param values (In) Matrix values.
     \param numEntries (In) Number of matrix coefficients being inserted.
  */
  int insertMatrixValues(int ownedPhysicsID, int ownedNode, int boxPhysicsID, int * boxNodeList, double * values, int numEntries);

  //! Insert matrix coefficients for a given row, where columns are from different physics types at the same node.
  /*! Insert values into matrix for a given row, where columns are from different physics types at the same node.
     \param ownedPhysicsID (In) 
                           This should be between 0 and one less than the number of physics variables tracked at a node.  
			   node and rowPhysicsID together specify which row of the matrix is being updated.
     \param ownedNode (In) Current owned node ID.  This is the local ID based on the set of owned nodes for this processor.
                      This should be between 0 and one less than the number of nodes owned by this processor, independent of
		      the number of physics types being computed for the given problem.
     \param boxPhysicsIDList (In) The list of indices for the types of unknowns to use for the column indices.  
                           This should be an array of length numEntries where each value is
			   between 0 and one less than the number of physics variables tracked at a node.  
			   boxPhysicsIDList and boxNode together specify which column of the matrix is being updated.
     \param boxNode (In) Nodal index corresponding to matrix values. This is a local ID in the box coordinate system.
     \param values (In) Matrix values.
     \param numEntries (In) Number of matrix coefficients being inserted.
  */
  int insertMatrixValues(int ownedPhysicsID, int ownedNode, int * boxPhysicsIDList, int boxNode, double * values, int numEntries);

  //! Method that must be called each time matrix value insertion is complete (usually once per nonlinear iteration).
  int finalizeProblemValues();

  //! Declare a block matrix to be read-only or read-write.
  /*! Change the read-only mode of a matrix block i,j.  Initially all blocks are read-write (not read-only).
     \param rowPhysicsID (In) ith physics block (row).
     \param colPhysicsID (In) jth physics block (column).
     \param readOnly (In) Bool:  If true, then the i,j block of the matrix will not be further modified.  
     Thus the method initializeProblemValues() will not zero out this block of the matrix.  Calling this method
     with readOnly set to false will allow initializeProblemValues() and insertMatrixValues() to change the contents
     of the i,j block.
  */
  int setBlockMatrixReadOnly(int rowPhysicsID, int colPhysicsID, bool readOnly);

  //@}

  //@{ \name Single-call modifier and accessor methods

  //! Set all left hand side (initial guess) vectors at once.
  /*! Allows the definition of initial guess values in a single call.
    \param b (In) An array of pointers of length numPhysicsTypes, where each array x[i] is of length numOwnedNodes.
  */
  int setRhs(const double ** b);

  //! Set all left hand side (solution) vectors at once.
  /*! Allows definition of all solution values in a single call.
    \param x (In) An array of pointers of length numPhysicsTypes, where each array x[i] is of length numBoxNodes.
  */
  int setLhs(double ** x) const;

  //! Get all right hand side vectors at once.
  /*! Allows the definition of initial guess values in a single call.
    \param b (Out) An array of pointers of length numPhysicsTypes, where each array x[i] is of length numOwnedNodes.
  */
  int getRhs(double ** b) const;

  //! Get all left hand side (solution) vectors at once.
  /*! Allows access to all solution values in a single call.
    \param x (Out) An array of pointers of length numPhysicsTypes, where each array x[i] is of length numBoxNodes.
  */
  int getLhs(double ** x) const;
  //@}

  //@{ \name Solver-related methods

  //! Setup up the solver for solving the current linear problem.
  int setupSolver();

  //! Solve the current linear problem.
  int solve();
  //@}

  //@{ \name Miscellaneous support methods (used by the application, or by LOCA, or both)

  //! Apply global linear operator for all physics types, b = Ax
  /*! Computes b = Ax where A is the global linear operator and each x[i] and b[i] is associated with the ith physics block.
    \param x (In) Array of pointers such that x[i] is an array of doubles of length numBoxNodes.
    \param b (Out) Array of pointers such that b[i] is an array of doubles of length numOwnedNodes.

    \warning Note that x[i] is of length numBoxNodes and b[i] is of length numOwnedNodes.
  */
  int applyMatrix(const double** x, double** b) const;
  
  //! Fill the arrays xBox[i] with xOwned[i] for all physics types i, i.e., fill in ghost values on each processor.
  /*! Fills in xBox[i] with values from xOwned[i] such that ghost values from other processors are updated.
    \param xOwned (In) Array of pointers such that xOwned[i] is an array of doubles of length numOwnedNodes.
    \param xBox (Out) Array of pointers such that xBox[i] is an array of doubles of length numBoxNodes.

    \warning Note that xOwned[i] is of length numOwnedNodes and xBox[i] is of length numBoxNodes.
  */
  int importR2C(const double** xOwned, double** xBox) const;

  //! Fill the array aBox with aOwned, i.e., fill in ghost values on each processor.
  //! Fill the arrays aBox with aOwned for a single physics types.
  /*! Fills in aBox with values from aOwned such that ghost values from other processors are updated.
    \param aOwned (In) An array of doubles of length numOwnedNodes.
    \param aBox (Out) An array of doubles of length numBoxNodes.

    \warning Note that aOwned is of length numOwnedNodes and aBox is of length numBoxNodes.
  */
  int importR2C(const double* aOwned, double* aBox) const;

  //! Fill the array aOwned with aBox, i.e., filter out ghost values on each processor.
  //! Fill the arrays aOwned with aBox for a single physics types.
  /*! Fills in aOwned with values from aBox such that ghost values from other processors are removed.
    \param aBox (In) An array of doubles of length numBoxNodes.
    \param aOwned (Out) An array of doubles of length numOwnedNodes.

    \warning Note that aOwned is of length numOwnedNodes and aBox is of length numBoxNodes.
  */
  int exportC2R(const double* aBox, double* aOwned) const;

  //@}

private:

  inline int ownedToSolverGID(int ownedPhysicsID, int ownedNode) { 
    if (groupByPhysics_) 
      return(ownedPhysicsID*numOwnedNodes + ownedMap.GID(ownedNode));
    else
      return(ownedPhysicsID + numOwnedNodes*ownedMap.GID(ownedNode);
  }
	     
  inline int ownedToSolverLID(int ownedPhysicsID, int ownedNode) { 
    if (groupByPhysics_) 
      return(ownedPhysicsID*numOwnedNodes + ownedNode);
    else
      return(ownedPhysicsID + numOwnedNodes*ownedNode);
  }
	     
  inline int boxToSolverGID(int boxPhysicsID, int boxNode) { 
    if (groupByPhysics_) 
      return(boxPhysicsID*numBoxNodes + boxMap.GID(boxNode));
    else
      return(boxPhysicsID + numBoxNodes*boxMap.GID(boxNode));
  }
	     
  int numUnknownsPerNode_;
  int numOwnedNodes_;
  int numBoxNodes_;
  int numMatrixBlocks_;
  Epetra_MpiComm comm_;
  Epetra_RowMatrix *** blockMatrix_;
  bool ** blockMatrixReadOnly_;
  bool ** blockMatrixIsVbr_;
  Epetra_CrsGraph *** blockGraph_;
  Epetra_Vector ** blockLhs_;
  Epetra_Vector ** blockRhs_;
  Epetra_BlockMap ** rowMaps_;
  Epetra_BlockMap ** colMaps_;
  Teuchos::RefCountPtr<Epetra_Map> ownedMap_;
  Teuchos::RefCountPtr<Epetra_Map> boxMap_;
  Teuchos::RefCountPtr<Epetra_Import> ownedToBoxImporter_;
  Teuchos::RefCountPtr<Epetra_CrsMatrix> globalMatrix_;
  Teuchos::RefCountPtr<Epetra_Vector> globalRhs_;
  Teuchos::RefCountPtr<Epetra_Vector> globalLhs_;
  Teuchos::RefCountPtr<AztecOO> solver_;
  bool isBlockStructureSet_;
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool groupByPhysics_;



};

#endif /* DFT_SOLVERMANAGER_HPP */
