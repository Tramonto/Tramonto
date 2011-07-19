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

class Epetra_Vector;
class Epetra_IntVector;
class Epetra_Map;
class Epetra_Import;
class Epetra_CrsMatrix;
class Epetra_LinearProblem;
class Epetra_Operator;
class AztecOO;
#include "Epetra_MpiComm.h"
#include "Teuchos_RefCountPtr.hpp"
#include <map>
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Teuchos_ParameterList.hpp"
#include "dft_direct_solver_const.h"

#include "Epetra_Vector.h"
#include "Epetra_Import.h"

//! dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>:  Solver manager class for Tramonto using Trilinos.
/*! The dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node> class supports solver capabilities for Tramonto.

*/

template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node=Kokkos::DefaultNode::DefaultNodeType>
class dft_BasicLinProbMgr {

  public:
  TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node);

  //@{ \name Constructors/destructor.

  //! dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node> Constructor.
  /*! Initialize a linear problem manager for Tramonto
     \param numUnknownsPerNode (In) The number of unknowns tracked per node of the mesh.
     \param solverOptions (In) An array of ints defined in dft_solver_defs.h containing information to guide and report solver status.
     \param solverParams (In) An array of Scalars defined in dft_solver_defs.h containing information to guide and report solver status.
     \param comm (In) Teuchos communicator that should be used by the solver.
  */
  dft_BasicLinProbMgr
  (MPI_Comm ecomm, size_t numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm);

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
  (int * GIDs, const ArrayView<const GlobalOrdinal>& GIDs2, LocalOrdinal nx = 0, LocalOrdinal ny = 1, LocalOrdinal nz = 1);

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
  (int * GIDs, const ArrayView<const GlobalOrdinal>& GIDs2, LocalOrdinal nx = 0, LocalOrdinal ny = 1, LocalOrdinal nz = 1);

  //! Define the nodes on this processor that will be mesh-coarsened, must be nodes set as part of setNodalRowMap().
  /*! Define the list of global node IDs that will be algebraically defined as an average of neighbors.  The exact definition
      of the averaging formula is not important.
     \param GIDs (In) List of global IDs that are flagged for coarsening.

     \warning If one processor calls this function, then ALL must call it.

  */
  virtual void
  setCoarsenedNodesList
  (const ArrayView<const GlobalOrdinal> &GIDs2);

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
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, const ArrayView<const LocalOrdinal> &boxNodeList, const ArrayView<const Scalar>& values);

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
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, const ArrayView<const LocalOrdinal> &boxPhysicsIDList, LocalOrdinal boxNode, const ArrayView<const Scalar> &values);

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
  (double ** x) const;
  //@}


  //@{ \name Solver-related methods

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
  (const double** xOwned, double** xBox, const ArrayRCP<const ArrayRCP<const Scalar> >& xOwned2) const;

  //! Fill the array aBox with aOwned, i.e., fill in ghost values on each processor.
  //! Fill the arrays aBox with aOwned for a single physics types.
  /*! Fills in aBox with values from aOwned such that ghost values from other processors are updated.
    \param aOwned (In) An ArrayRCP of Scalars of length numOwnedNodes.
    \param aBox (Out) An ArrayRCP of Scalars of length numBoxNodes.

    \warning Note that aOwned is of length numOwnedNodes and aBox is of length numBoxNodes.
  */
  virtual ArrayRCP<Scalar>
  importR2C
  (const double* aOwned, double* aBox, const ArrayRCP<const Scalar> &aOwned2) const;

  //! Fill the array aOwned with aBox, i.e., filter out ghost values on each processor.
  //! Fill the arrays aOwned with aBox for a single physics types.
  /*! Fills in aOwned with values from aBox such that ghost values from other processors are removed.
    \param aBox (In) An ArrayRCP of Scalars of length numBoxNodes.
    \param aOwned (Out) An ArrayRCP of Scalars of length numOwnedNodes.

    \warning Note that aOwned is of length numOwnedNodes and aBox is of length numBoxNodes.
  */
  virtual ArrayRCP<Scalar>
  exportC2R
  (const ArrayRCP<const Scalar>& aBox2) const;
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
  RCP<MAT> globalMatrix_;
  RCP<VEC> globalRhs_;
  RCP<VEC> globalLhs_;
  RCP<ParameterList> parameterList_;
  GlobalOrdinal curRow_;
  std::map<GlobalOrdinal, Scalar> curRowValues_;
  Array<GlobalOrdinal> indices_;
  Array<Scalar> values_;

  // Epetra variables
  Epetra_MpiComm commEP_;
  Teuchos::RefCountPtr<Epetra_Map> ownedMapEP_;
  Teuchos::RefCountPtr<Epetra_Map> boxMapEP_;
  Teuchos::RefCountPtr<Epetra_Import> ownedToBoxImporterEP_;

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
  RCP<PRECOND> preconditioner_;
#endif

};

#endif

//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
dft_BasicLinProbMgr
(MPI_Comm ecomm, size_t numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm)
  : commEP_(Epetra_MpiComm(ecomm)),
    isBlockStructureSet_(false),
    isGraphStructureSet_(false),
    isLinearProblemSet_(false),
    groupByPhysics_(true),
    firstTime_(true),
    numUnknownsPerNode_(numUnknownsPerNode),
    parameterList_(parameterList),
    numOwnedNodes_(0),
    numBoxNodes_(0),
    numGlobalNodes_(0),
    numGlobalBoxNodes_(0),
    numCoarsenedNodes_(0),
    numGlobalCoarsenedNodes_(0),
    comm_(comm),
    curRow_(-1) {
    return;
    }
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
~dft_BasicLinProbMgr()
{
  return;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalRowMap
(int *GIDs, const ArrayView<const GlobalOrdinal>& GIDs2, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz)
{
  if (numGlobalNodes_!=0)
  {
    return; // Already been here
  }
  numOwnedNodes_ = GIDs2.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
				  &numOwnedNodes_, &numGlobalNodes_);

  ownedMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs2, 0, comm_));

  // Epetra code
  ownedMapEP_ = Teuchos::rcp(new Epetra_Map(-1, numOwnedNodes_, GIDs, 0, commEP_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setNodalColMap
(int *GIDs, const ArrayView<const GlobalOrdinal> &GIDs2, LocalOrdinal nx, LocalOrdinal ny, LocalOrdinal nz)
{
  if (numGlobalBoxNodes_!=0)
  {
    return; // Already been here
  }

  numBoxNodes_ = GIDs2.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
			       &numBoxNodes_, &numGlobalBoxNodes_);

  boxMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs2, 0, comm_));

  // Epetra code
  boxMapEP_ = Teuchos::rcp(new Epetra_Map(-1, numBoxNodes_, GIDs, 0, commEP_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setCoarsenedNodesList(const ArrayView<const GlobalOrdinal> &GIDs)
{
  if (numGlobalCoarsenedNodes_!=0)
  {
    return; // Already been here
  }

  numCoarsenedNodes_ = GIDs.size();
  Teuchos::reduceAll<int, size_t>(*comm_, Teuchos::REDUCE_SUM, 1,
			       &numCoarsenedNodes_, &numGlobalCoarsenedNodes_);

  coarsenedNodesMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), GIDs, 0, comm_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeBlockStructure
()
{
  TEST_FOR_EXCEPTION(isBlockStructureSet_, std::runtime_error, "Already set block structure.\n");

  const size_t numUnks2 = numOwnedNodes_*numUnknownsPerNode_;
  Array<LocalOrdinal> globalGIDList2(numUnks2);

  // Physics ordering for Basic Linear Problem is natural ordering:
  physicsOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++)
  {
    physicsOrdering_[i] = i;
  }

  // create inverse mapping of where each physics unknown is ordered for the solver
  solverOrdering_.resize(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<physicsOrdering_.size(); i++)
  {
    solverOrdering_[physicsOrdering_[i]]=i;
  }

  // Sanity check of physics ordering
  checkPhysicsOrdering();

  ArrayView<const GlobalOrdinal> GIDs2 = ownedMap_->getNodeElementList();
  LocalOrdinal k2=0;
  if (groupByPhysics_)
  {
    for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++)
    {
      LocalOrdinal ii=physicsOrdering_[i];
      for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
      {
	      globalGIDList2[k2++] = ii*numGlobalNodes_ + GIDs2[j];
      }
    }
  }
  else
  {
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++)
    {
      for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++)
      {
	LocalOrdinal ii=physicsOrdering_[i];
	globalGIDList2[k2++] = ii + GIDs2[j]*numUnknownsPerNode_;
      }
    }
  }

  globalRowMap_ = rcp(new MAP(Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), globalGIDList2, 0, comm_));
  globalMatrix_ = rcp(new MAT(globalRowMap_, 0));
  globalMatrix_->setObjectLabel("BasicLinProbMgr::globalMatrix");
  globalRhs_ = rcp(new VEC(globalRowMap_));
  globalLhs_ = rcp(new VEC(globalRowMap_));
  ownedToBoxImporter_ = rcp(new Tpetra::Import<LocalOrdinal, GlobalOrdinal, Node>(ownedMap_, boxMap_));

  isBlockStructureSet_ = true;
  isGraphStructureSet_ = true;

  //// Epetra Code ////
  ownedToBoxImporterEP_ = Teuchos::rcp(new Epetra_Import(*boxMapEP_, *ownedMapEP_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
initializeProblemValues
()
{
  TEST_FOR_EXCEPTION(!isBlockStructureSet_, std::logic_error,
		     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");
  TEST_FOR_EXCEPTION(!isGraphStructureSet_, std::logic_error,
		     "Linear problem structure must be completely set up.  This requires a sequence of calls, ending with finalizeBlockStructure");

  isLinearProblemSet_ = false; // We are reinitializing the linear problem

 // AGS: I found that we needed to initialize the matrix even the
 // first time a matrix was filled, because matrix entries are being put
 // in on residual-only fills, which can occur before matrix fills.

  if (!firstTime_)
  {
    globalMatrix_->resumeFill();
    globalMatrix_->setAllToScalar(0.0);
    globalMatrix_->fillComplete();
    globalRhs_->putScalar(0.0);
    globalLhs_->putScalar(0.0);
  }
 }
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRhsValue
(GlobalOrdinal ownedPhysicsID, GlobalOrdinal ownedNode, Scalar value)
{
  LocalOrdinal rhsLID = ownedToSolverLID(ownedPhysicsID, ownedNode); // Get solver LID
  globalRhs_->sumIntoLocalValue(rhsLID, value);
  //cout << std::setprecision(2);
  //cout << "b[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode<<"] = " << value << endl;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, Scalar value)
{
  GlobalOrdinal rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = boxToSolverGID(boxPhysicsID, boxNode);

  //cout << std::setprecision(2);
  //cout << "A[ownedPhysicsID="<<ownedPhysicsID<<"][ownedNode="<<ownedNode
  //     << "][boxPhysicsID="  <<boxPhysicsID  <<"][boxNode="  <<boxNode
  //     << "][rowGID="        <<rowGID        <<"][colGID="   <<colGID
  //     << "] = " << value << endl;
  if (firstTime_) {
    if (rowGID!=curRow_) {
      insertRow();  // Dump the current contents of curRowValues_ into matrix and clear map
      curRow_=rowGID;
    }
    curRowValues_[colGID] += value;
  }
  else {
    Array<Scalar> vals(1);
    vals[0] = value;
    Array<GlobalOrdinal> globalCols(1);
    globalCols[0] = colGID;
    globalMatrix_->sumIntoGlobalValues(rowGID, globalCols, vals);
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertRow
()
{
  if (curRowValues_.empty()) return;

  typename std::map<LocalOrdinal, Scalar>::iterator pos2;
  for (pos2 = curRowValues_.begin(); pos2 != curRowValues_.end(); ++pos2) {
    indices_.append(pos2->first);
    values_.append(pos2->second);
  }

  globalMatrix_->insertGlobalValues(curRow_, indices_, values_);

  indices_.clear();
  values_.clear();
  curRowValues_.clear();
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Scalar
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getMatrixValue
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode)
{
  TEST_FOR_EXCEPTION(globalMatrix_.get()==0, std::runtime_error, "Global Matrix is not constructed, must set debug flag to enable this feature.\n");

  GlobalOrdinal rowGID = ownedToSolverGID(ownedPhysicsID, ownedNode); // Get solver Row GID
  GlobalOrdinal colGID = boxToSolverGID(boxPhysicsID, boxNode);
  size_t numEntries;
  ArrayView<const GlobalOrdinal> indices;
  ArrayView<const Scalar> values;
  if (globalMatrix_->isGloballyIndexed())
  {
    globalMatrix_->getGlobalRowView(rowGID, indices, values); // get view of current row
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; i++)
    {
      if (colGID==indices[i])
      {
	return(values[i]);
      }
    }
  }
  else
  {
    rowGID = globalMatrix_->getRowMap()->getLocalElement(rowGID); // get local row ID
    colGID = globalMatrix_->getColMap()->getLocalElement(colGID); // get local column ID
    globalMatrix_->getLocalRowView(rowGID, indices, values); // get view of current row
    numEntries = indices.size();
    for (LocalOrdinal i=0; i<numEntries; i++)
    {
      if (colGID==indices[i])
      {
	return(values[i]);
      }
    }
  }

  return(0.0);
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID,
 const ArrayView<const LocalOrdinal>& boxNodeList, const ArrayView<const Scalar>& values)
{
  for (LocalOrdinal i=0; i<boxNodeList.size(); i++)
  {
    insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsID, boxNodeList[i], values[i]);
  } //end for
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
insertMatrixValues
(LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, const ArrayView<const LocalOrdinal> &boxPhysicsIDList,
 LocalOrdinal boxNode, const ArrayView<const Scalar> &values)
{
  for (LocalOrdinal i=0; i<boxPhysicsIDList.size(); i++)
  {
    insertMatrixValue(ownedPhysicsID, ownedNode, boxPhysicsIDList[i], boxNode, values[i]);
  } //end for
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
finalizeProblemValues
()
{
  if (isLinearProblemSet_)
    return; // nothing to do

  if (firstTime_) {
    insertRow();

    if(!globalMatrix_->isFillComplete()){
      globalMatrix_->fillComplete();
    }
  }

  isLinearProblemSet_ = true;
  firstTime_ = false;

  //writeMatrix("basica.dat", "", "");
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setRhs
(const ArrayView<const ArrayView<const Scalar> >& b)
{
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      globalRhs_->replaceLocalValue(ownedToSolverLID(i,j), b[i][j]);
    }
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setLhs
(const ArrayView<const ArrayView<const Scalar> > &x) const
{
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    ArrayRCP<const Scalar> xtmp = exportC2R(Teuchos::arcpFromArrayView<const Scalar>(x[i])); // Use simple import
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      globalLhs_->replaceLocalValue(ownedToSolverLID(i,j), xtmp[j]);
    }
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getLhs
(double ** x) const
{
  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalLhs_->get1dViewNonConst();
  Epetra_SerialDenseVector xtmp(numOwnedNodes_); // Temp vector to hold local x values

  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    ArrayRCP<Scalar> temp(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      xtmp[j] = tmp[ownedToSolverLID(i,j)]; // Epetra code
      temp[j] = tmp[ownedToSolverLID(i,j)];
    }
    ArrayOfPtrs[i] = importR2C(xtmp.Values(),x[i],temp.getConst()); // Use simple import (uses view so doesn't work)
  }
  return ArrayOfPtrs;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
getRhs
() const
{
  ArrayRCP<ArrayRCP<Scalar> > ArrayOfPtrs = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  ArrayRCP<Scalar> tmp = globalRhs_->get1dViewNonConst();
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    ArrayOfPtrs[i] = Teuchos::arcp<Scalar>(numOwnedNodes_);
    for (LocalOrdinal j=0; j<numOwnedNodes_; j++) {
      ArrayOfPtrs[i][j] = tmp[ownedToSolverLID(i,j)];
    }
  }
  return ArrayOfPtrs;
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
setupSolver
()
{
  TEST_FOR_EXCEPTION(!isLinearProblemSet_, std::logic_error,
		     "Linear problem must be completely set up.  This requires a sequence of calls, ending with finalizeProblemValues");

#ifdef SUPPORTS_STRATIMIKOS
  thyraRhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalRhs_);
  thyraLhs_ = createVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalLhs_);
  thyraOp_ = createLinearOp<Scalar, LocalOrdinal, GlobalOrdinal, Node>(globalMatrix_);

  solver_ = rcp(new DefaultLinearSolverBuilder("./dft_input.xml"));
  RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
  solver_->readParameters(out.get());

  lowsFactory_ = solver_->createLinearSolveStrategy("");
  lows_ = linearOpWithSolve<Scalar>(*lowsFactory_, thyraOp_);
#else
  //printf("Node num rows = %d\n", globalMatrix_->getNodeNumRows());
  //printf("Node num cols = %d\n", globalMatrix_->getNodeNumCols());
  problem_ = rcp(new LinPROB(globalMatrix_, globalLhs_, globalRhs_));
  // No preconditioner for now
  //  preconditioner_ = Ifpack2::Factory::create<MAT>("ILUT", globalMatrix_);
  //  preconditioner_->setParameters(*parameterList_);
  //  preconditioner_->initialize();
  //  preconditioner_->compute();
  //  problem_->setLeftPrec(preconditioner_);
  TEST_FOR_EXCEPT(problem_->setProblem() == false);
  solver_ = rcp(new Belos::BlockGmresSolMgr<Scalar, MV, OP>(problem_, parameterList_));
#endif
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
solve
()
{
#ifdef KDEBUG
  printf("\n\n\n\ndft_BasicLinProbMgr::solve()\n\n\n\n");
#endif

#ifdef SUPPORTS_STRATIMIKOS
  SolveStatus<Scalar> status = lows_->solve(Thyra::NOTRANS, *thyraRhs_, thyraLhs_.ptr());
#else
  ReturnType ret = solver_->solve();
#endif

  bool writeMatrixNow = false;
  if (writeMatrixNow)
  {
    writeMatrix("A_ref.dat", "GlobalMatrix", "GlobalMatrix");
    writeLhs("x_ref.dat");
    writeRhs("b_ref.dat");
    writePermutation("p_ref.dat");
  }
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
applyMatrix
(const ArrayView<const ArrayView<const Scalar> >& x) const
{
  setLhs(x);
  globalMatrix_->apply(*globalLhs_, *globalRhs_);
  return(getRhs());
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<ArrayRCP<Scalar> >
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const double** xOwned, double** xBox, const ArrayRCP<const ArrayRCP<const Scalar> >& xOwned2) const
{
  ArrayRCP<ArrayRCP<Scalar> > my_xBox = Teuchos::arcp<ArrayRCP<Scalar> >(numUnknownsPerNode_);
  for (LocalOrdinal i=0; i<numUnknownsPerNode_; i++) {
    my_xBox[i] = importR2C(xOwned[i], xBox[i], xOwned2[i]);
  }

  return(my_xBox);
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
importR2C
(const double* aOwned, double* aBox, const ArrayRCP<const Scalar> &aOwned2) const
{
  RCP<VEC> owned2 = rcp(new VEC(ownedMap_, aOwned2()));
  RCP<VEC>  box2 = rcp(new VEC(boxMap_));
  box2->doImport(*owned2, *ownedToBoxImporter_, Tpetra::INSERT);

  // Epetra code
  Epetra_Vector owned(View, *ownedMapEP_, (double *) aOwned);
  Epetra_Vector box(View, *boxMapEP_, aBox);
  box.Import(owned, *ownedToBoxImporterEP_, Insert);

  return (box2->get1dViewNonConst());
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
ArrayRCP<Scalar>
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
exportC2R
(const ArrayRCP<const Scalar>& aBox) const
{
  RCP<VEC> owned2 =  rcp(new VEC(ownedMap_));
  RCP<VEC> box = rcp(new VEC(boxMap_, aBox()));
  owned2->doExport(*box, *ownedToBoxImporter_, Tpetra::INSERT); // Use importer, but zero out off-processor contributions.

  return(owned2->get1dViewNonConst());
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeMatrix
(const char * filename, const char * matrixName, const char * matrixDescription) const  {
  //int dft_BasicLinProbMgr::writeMatrix(const char * filename, const char * matrixName, const char * matrixDescription) const  {
  return;
      //(EpetraExt::RowMatrixToMatrixMarketFile(filename, *globalMatrix_, matrixName, matrixDescription));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeLhs
(const char * filename) const  {
  //int dft_BasicLinProbMgr::writeLhs(const char * filename) const  {
  return;
      //(EpetraExt::MultiVectorToMatlabFile(filename, *globalLhs_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writeRhs
(const char * filename) const  {
  //int dft_BasicLinProbMgr::writeRhs(const char * filename) const  {
  return;
      //(EpetraExt::MultiVectorToMatlabFile(filename, *globalRhs_));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
writePermutation
(const char * filename) const  {
  //int dft_BasicLinProbMgr::writePermutation(const char * filename) const  {
  return;
    //(EpetraExt::BlockMapToMatrixMarketFile(filename, *globalRowMap_, " ", " ", false));
}
//=============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void
dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
checkPhysicsOrdering()
const  {
  TEST_FOR_EXCEPTION(physicsOrdering_.size()==0, std::runtime_error, "No unknowns are registered with this problem manager.\n");

  size_t numUnks = physicsOrdering_.size();
  Array<Scalar> tmp(numUnks);
  for (LocalOrdinal i=0; i<numUnks; i++)
  {
    LocalOrdinal curID = physicsOrdering_[i];
    TEST_FOR_EXCEPTION(curID <0, std::runtime_error, "Invalid unknown number " << curID << " is less than 0.\n");
    TEST_FOR_EXCEPTION(curID>=numUnks, std::runtime_error, "Invalid unknown number " << curID << " is greater than or equal to the number of unknowns (" << numUnks << ").\n");
      tmp[curID] = tmp[curID]+1;
      // Increment counter for this ID (at the end each ID should appear exactly one time).
  }

  for (LocalOrdinal i=0; i<numUnks; i++)
  {
    TEST_FOR_EXCEPTION(tmp[i]==0, std::runtime_error, "Unknown number " << i << " is not present and should be.\n");
    TEST_FOR_EXCEPTION(tmp[i]>1, std::runtime_error, "Unknown number " << i << " is present " << tmp[i] << " times and should be present only once.\n");
  }

}
template class dft_BasicLinProbMgr<double, int, int>;
