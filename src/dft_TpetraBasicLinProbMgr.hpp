/*@HEADER
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
*/
#ifndef DFT_BASICLINPROBMGR_HPP
#define DFT_BASICLINPROBMGR_HPP

#include "Tpetra_Headers.hpp"
#include "dft_direct_solver_const.h"

//! dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>:  Solver manager class for Tramonto using Trilinos.
/*! The dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node> class supports solver capabilities for Tramonto.
*/

template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node=Kokkos::DefaultNode::DefaultNodeType>
class dft_BasicLinProbMgr {

  public:
  TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node)
  TYPEDEF_MIXED(Scalar, LocalOrdinal, GlobalOrdinal, Node)

  //@{ \name Constructors/destructor.

  //! dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node> Constructor.
  /*! Initialize a linear problem manager for Tramonto
     \param numUnknownsPerNode (In) The number of unknowns tracked per node of the mesh.
     \param solverOptions (In) An array of ints defined in dft_solver_defs.h containing information to guide and report solver status.
     \param solverParams (In) An array of Scalars defined in dft_solver_defs.h containing information to guide and report solver status.
     \param comm (In) Teuchos communicator that should be used by the solver.
  */
  dft_BasicLinProbMgr
  (size_t numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm);

  //! dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node> Destructor.
  /*! Completely deletes a dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node> object.
  */
  virtual ~dft_BasicLinProbMgr();
  //@}


  //@{ \name Block structure setup methods

  //! Define global to local index row mappings
  /*! Define the list of global node IDs owned by this processor.  This mapping is used for identifying
     the global row indexing given the local row IDs and the physics ID.
     \param numOwnedNodes (In) Number of owned nodes.
     \param GIDs (In) List of global IDs in the order of local indexing of the owned nodes.
     \param nx (In) Number of nodes in X direction on this processor.  If set to 0 (default value), nx, ny and nz will be ignored.  For a 1D problem, this value should be set to numGIDs.
     \param ny (In) Number of nodes in Y direction on this processor, defaults to 1.
     \param nz (In) Number of nodes in Z direction on this processor, defaults to 1.

  */
  virtual void
  setNodalRowMap
  (const ArrayView<const GlobalOrdinal>& GIDs, LocalOrdinal nx = 0, LocalOrdinal ny = 1, LocalOrdinal nz = 1);

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
  virtual void
  setNodalColMap
  (const ArrayView<const GlobalOrdinal>& GIDs, LocalOrdinal nx = 0, LocalOrdinal ny = 1, LocalOrdinal nz = 1);

  //! Define the nodes on this processor that will be mesh-coarsened, must be nodes set as part of setNodalRowMap().
  /*! Define the list of global node IDs that will be algebraically defined as an average of neighbors.  The exact definition
      of the averaging formula is not important.
     \param GIDs (In) List of global IDs that are flagged for coarsening.

     \warning If one processor calls this function, then ALL must call it.

  */
  virtual void
  setCoarsenedNodesList
  (const ArrayView<const GlobalOrdinal> &GIDs);

  //! Method that must be called once, when all row and column maps are set.
  /*! This method constructs all of the CrsGraph objects and the lhs and rhs vectors. */
  virtual void
  finalizeBlockStructure
  ();
  //@}


  //@{ \name Matrix, lhs and rhs value setup methods for piece-wise construction.

  //! Method that must be called each time \e prior to starting matrix, lhs and rhs value insertion (usually once per nonlinear iteration).
  /*! This method zeros out the matrix, lhs and rhs values. */
  virtual void
  initializeProblemValues
  ();

  //! Insert rhs value based on ownedNode and ownedPhysicsID.
  /*! Insert rhs value into entry based on owned node and physicsID.
     \param ownedPhysicsID (In) The index for the type of unknown.
			   This should be between 0 and one less than the number of physics variables tracked at a node.
     \param ownedNode (In) Current owned node ID.  This is the local ID based on the set of owned nodes for this processor.
		      This should be between 0 and one less than the number of nodes owned by this processor, independent of
		      the number of physics types being computed for the given problem.
     \param value (In) Rhs value.
  */
  virtual void
  insertRhsValue
  (GlobalOrdinal ownedPhysicsID, GlobalOrdinal ownedNode, Scalar value);

  //! Insert single matrix coefficient into system.
  /*! Insert rhs value into entry based on owned node and physicsID.
     \param ownedPhysicsID (In) The index for the type of unknown.
			   This should be between 0 and one less than the number of physics variables tracked at a node.
     \param ownedNode (In) Current owned node ID.  This is the local ID based on the set of owned nodes for this processor.
		      This should be between 0 and one less than the number of nodes owned by this processor, independent of
		      the number of physics types being computed for the given problem.
     \param value (In) Rhs value.
  */
  virtual void
  insertMatrixValue
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, Scalar value);

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
  */
  virtual void
  insertMatrixValues
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID,
   const ArrayView<const LocalOrdinal> &boxNodeList, const ArrayView<const Scalar>& values);

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
  */
  virtual void
  insertMatrixValues
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, const ArrayView<const LocalOrdinal> &boxPhysicsIDList,
   LocalOrdinal boxNode, const ArrayView<const Scalar> &values);

  //! Method that must be called each time matrix value insertion is complete (usually once per nonlinear iteration).
  virtual void
  finalizeProblemValues
  ();
  //@}


  //@{ \name Single-call modifier and accessor methods

  //! Get a matrix entry
  /*! Returns the matrix value for the given physics and node values, if the value is not present then a zero value is returned.
     \param ownedPhysicsID (In) The physics ID for the matrix row being queried.
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
     \return The requested matrix value is returned; if the value is not present, the return value will be zero.
  */
  virtual Scalar
  getMatrixValue
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode);

  //! Set all left hand side (initial guess) vectors at once.
  /*! Allows the definition of initial guess values in a single call.
    \param b (In) An ArrayView of ArrayViews of length numPhysicsTypes, where each array x[i] is of length numOwnedNodes.
  */
  virtual void
  setRhs
  (const ArrayView<const ArrayView<const Scalar> >& b);

  //! Set all left hand side (solution) vectors at once.
  /*! Allows definition of all solution values in a single call.
    \param x (In) An ArrayView of ArrayViews of length numPhysicsTypes, where each array x[i] is of length numBoxNodes.
  */
  virtual void
  setLhs
  (const ArrayView<const ArrayView<const Scalar> > &x) const;

  //! Get all right hand side vectors at once.
  /*! Allows the definition of initial guess values in a single call.
    \param b (Out) An ArrayRCP of ArrayRCPs of length numPhysicsTypes, where each array x[i] is of length numOwnedNodes.
  */
  virtual ArrayRCP<ArrayRCP<Scalar> >
  getRhs
  () const;

  //! Get all left hand side (solution) vectors at once.
  /*! Allows access to all solution values in a single call.
    \param x (Out) An ArrayRCP of ArrayRCPs of length numPhysicsTypes, where each array x[i] is of length numBoxNodes.
  */
  virtual ArrayRCP<ArrayRCP<Scalar> >
  getLhs
  () const;
  //@}


  //@{ \name Solver-related methods

  //! Get constants for this machine and problem.
  virtual void
  setMachineParams
  ();

  //! Setup up the solver for solving the current linear problem.
  virtual void
  setupSolver
  ();

  //! Solve the current linear problem.
  virtual void
  solve
  ();
  //@}


  //@{ \name Output facilities

  //! Write matrix to specified filename using Matrix Market (i,j,value) format.
  /*! Write the matrix owned by the solver manager to the file called filename on PE 0 regardless of how many processors are
      involved in the computation.
      \param filename (In) Name of the file.  Any contents in this file will be erased.
      \param matrixName (In) Optional brief name for the matrix, will be inserted in the header of the output file, can be set to 0.
      \param matrixDescription (In) Optional longer description for the matrix, will be inserted in the header of the output file, can be set to 0.
  */
  virtual void
  writeMatrix
  (const char * filename, const char * matrixName, const char * matrixDescription) const;

  //! Write right hand side to specified filename in Matlab-compatible format.
  virtual void
  writeRhs
  (const char * filename) const;

  //! Write left hand side to specified filename in Matlab-compatible format.
  virtual void
  writeLhs
  (const char * filename) const;

  //! Write the permutation applied to the problem in Matlab-compatible format.
  virtual void
  writePermutation
  (const char * filename) const;
  //@}


  //@{ \name Miscellaneous support methods (used by the application, or by LOCA, or both)
  //! Apply global linear operator for all physics types, b = Ax
  /*! Computes b = Ax where A is the global linear operator and each x[i] and b[i] is associated with the ith physics block.
    \param x (In) ArrayView of ArrayViews such that x[i] is an array of Scalars of length numBoxNodes.
    \param b (Out) ArrayRCP of ArrayRCPs such that b[i] is an array of Scalars of length numOwnedNodes.

    \warning Note that x[i] is of length numBoxNodes and b[i] is of length numOwnedNodes.
  */
  virtual ArrayRCP<ArrayRCP<Scalar> >
  applyMatrix
  (const ArrayView<const ArrayView<const Scalar> >& x) const;

  //! Fill the arrays xBox[i] with xOwned[i] for all physics types i, i.e., fill in ghost values on each processor.
  /*! Fills in xBox[i] with values from xOwned[i] such that ghost values from other processors are updated.
    \param xOwned (In) ArrayRCP of ArrayRCPs such that xOwned[i] is an array of Scalars of length numOwnedNodes.
    \param xBox (Out) ArrayRCP of ArrayRCPs such that xBox[i] is an array of Scalars of length numBoxNodes.

    \warning Note that xOwned[i] is of length numOwnedNodes and xBox[i] is of length numBoxNodes.
  */
  virtual ArrayRCP<ArrayRCP<Scalar> >
  importR2C
  (const ArrayRCP<const ArrayRCP<const Scalar> >& xOwned) const;

  //! Fill the array aBox with aOwned, i.e., fill in ghost values on each processor.
  //! Fill the arrays aBox with aOwned for a single physics types.
  /*! Fills in aBox with values from aOwned such that ghost values from other processors are updated.
    \param aOwned (In) An ArrayRCP of Scalars of length numOwnedNodes.
    \param aBox (Out) An ArrayRCP of Scalars of length numBoxNodes.

    \warning Note that aOwned is of length numOwnedNodes and aBox is of length numBoxNodes.
  */
  virtual ArrayRCP<Scalar>
  importR2C
  (const ArrayRCP<const Scalar> &aOwned) const;

  //! Fill the array aOwned with aBox, i.e., filter out ghost values on each processor.
  //! Fill the arrays aOwned with aBox for a single physics types.
  /*! Fills in aOwned with values from aBox such that ghost values from other processors are removed.
    \param aBox (In) An ArrayRCP of Scalars of length numBoxNodes.
    \param aOwned (Out) An ArrayRCP of Scalars of length numOwnedNodes.

    \warning Note that aOwned is of length numOwnedNodes and aBox is of length numBoxNodes.
  */
  virtual ArrayRCP<Scalar>
  exportC2R
  (const ArrayRCP<const Scalar>& aBox) const;
  //@}


  virtual size_t
  getNumUnknownsPerNode
  ()
  {
    return numUnknownsPerNode_;
  }

  virtual size_t
  getNumOwnedNodes
  ()
  {
    return numOwnedNodes_;
  }

  virtual size_t
  getNumBoxNodes
  ()
  {
    return numBoxNodes_;
  }

protected:

  void
  checkPhysicsOrdering
  () const;

  inline GlobalOrdinal ownedToSolverGID(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode) const
  {
    if (groupByPhysics_)
    {
      return(ownedPhysicsID*numGlobalNodes_ + ownedMap_->getGlobalElement(ownedNode));
    }
    else
    {
      return(ownedPhysicsID + numUnknownsPerNode_*ownedMap_->getGlobalElement(ownedNode));
    }
  }

  inline
  LocalOrdinal
  ownedToSolverLID
  (GlobalOrdinal ownedPhysicsID, LocalOrdinal ownedNode) const
  {
    return(globalRowMap_->getLocalElement(ownedToSolverGID(ownedPhysicsID,ownedNode)));
  }

  inline int boxToSolverGID(int boxPhysicsID, int boxNode) const {
    if (groupByPhysics_)
      return(boxPhysicsID*numGlobalNodes_ + boxMap_->getGlobalElement(boxNode));
    else
      return(boxPhysicsID + numUnknownsPerNode_*boxMap_->getGlobalElement(boxNode));
  }

protected:
  void insertRow();
  bool isBlockStructureSet_;
  bool isGraphStructureSet_;
  bool isLinearProblemSet_;
  bool machineParamsSet_;
  bool groupByPhysics_;
  bool firstTime_;
  size_t numUnknownsPerNode_;
  size_t numOwnedNodes_;
  size_t numBoxNodes_;
  size_t numGlobalNodes_;
  size_t numGlobalBoxNodes_;
  size_t numCoarsenedNodes_;
  size_t numGlobalCoarsenedNodes_;
  RCP<const COMM> comm_;
  Array<GlobalOrdinal> physicsOrdering_;
  Array<GlobalOrdinal> solverOrdering_;
  RCP<const MAP> ownedMap_;
  RCP<const MAP> boxMap_;
  RCP<const MAP> coarsenedNodesMap_;
  RCP<VEC> ownedNodeIsCoarsened_;
  RCP<VEC> boxNodeIsCoarsened_;
  RCP<IMP> ownedToBoxImporter_;
  RCP<const MAP> globalRowMap_;
  RCP<MAT_P> globalMatrix_;
  RCP<OP> globalMatrixOperator_;
  RCP<SCALE_P> scalingMatrix_;
  RCP<VEC_P> rowScaleFactors_;
  RCP<VEC> rowScaleFactorsScalar_;
  RCP<VEC> globalRhs_;
  RCP<VEC> globalLhs_;
  RCP<ParameterList> parameterList_;
  RCP<ParameterList> tpetraParameterList_;
  GlobalOrdinal curRow_;
  Array<GlobalOrdinal> indices_;
  std::map<GlobalOrdinal, precScalar> curRowValues_;
  Array<precScalar> values_;
  GlobalOrdinal scaling_;
  GlobalOrdinal n_;
  Scalar eps_;
  halfScalar epsHalf_;
  Scalar anorm_;
  Scalar nae_;
  Scalar snae_;
  halfScalar naeHalf_;
  halfScalar snaeHalf_;

#ifdef SUPPORTS_STRATIMIKOS
  //Thyra Objects
  RCP<ThyraVEC> thyraRhs_;
  RCP<ThyraVEC> thyraLhs_;
  RCP<ThyraOP> thyraOp_;
  RCP<ThyraOP> precOp_;
  RCP<DefaultLinearSolverBuilder> solver_;

  RCP<ThyraLOWSFactory> lowsFactory_;
  RCP<ThyraLOWS> lows_;
  RCP<ThyraPRECFactory> precFactory_;
  RCP<ThyraPREC> prec_;
#else
  RCP<SolMGR> solver_;
  RCP<LinPROB> problem_;
  RCP<PRECOND_AS> preconditioner_;
  RCP<PRECOND_AS_OP> preconditionerOp_;
#endif

};

#endif
