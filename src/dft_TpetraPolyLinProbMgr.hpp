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

#ifndef DFT_POLYLINPROBMGR_HPP
#define DFT_POLYLINPROBMGR_HPP

#include "Tpetra_Headers.hpp"
#include "dft_PolyA11_Tpetra_Operator.hpp"
#include "dft_PolyA11_Coulomb_Tpetra_Operator.hpp"
#include "dft_PolyA22_Tpetra_Operator.hpp"
#include "dft_PolyA22_Coulomb_Tpetra_Operator.hpp"
#include "dft_Schur_Tpetra_Operator.hpp"
#include "dft_TpetraBasicLinProbMgr.hpp"

//====================================================

//! dft_PolyLinProbMgr:  Problem manager class for polymer problems.
/*! The dft_PolyLinProbMgr class supports polymer solver capabilities for Tramonto.

*/
template<class Scalar,class MatScalar=Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node=Kokkos::DefaultNode::DefaultNodeType>
class dft_PolyLinProbMgr:
  public virtual dft_BasicLinProbMgr<Scalar, MatScalar, LocalOrdinal, GlobalOrdinal, Node>
{
  public:
TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node)
TYPEDEF_MIXED(Scalar, LocalOrdinal, GlobalOrdinal, Node)

  typedef dft_BasicLinProbMgr<Scalar,MatScalar,LocalOrdinal,GlobalOrdinal,Node> BLPM;
  typedef dft_PolyA11_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> P11TO;
  typedef dft_PolyA11_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal> P11CO;
  typedef dft_PolyA22_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> P22TO;
  typedef dft_PolyA22_Coulomb_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal> P22CO;
  typedef dft_Schur_Tpetra_Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> ScTO;

  //@{ \name Constructors/destructor.
  //! dft_PolyLinProbMgr Constructor.
  /*! Initialize a linear problem manager for Tramonto
     \param numUnknownsPerNode (In) The number of unknowns tracked per node of the mesh.
     \param parameterList (In) A ParameterList containing information to guide and report solver status.
     \param comm (In) Teuchos communicator that should be used by the solver.
     \param debug (In) Turns debug mode on if set to true, false by default.
  */
  dft_PolyLinProbMgr(LocalOrdinal numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm, bool debug = false);

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
  void
  setGEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (gEquations_.size()!=0)
    {
      return; // Already been here
    } //end if
    gEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++)
    {
      gEquations_[i] = physicsIDs[i];
    } //end for
  } //end setGEquationIDs

  //! Define G inverse equation IDs
  /*! Define the list of physics IDs associated with the G inverse equations in descending order.
     \param numGInvEquations (In) Number of G inverse equations per node.
     \param physicsIDs (In) List of physics IDs associated with the G equations in ascending order.
  */
  void
  setGInvEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (gInvEquations_.size()!=0)
    {
      return; // Already been here
    } //end if
    gInvEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++)
    {
      gInvEquations_[i] = physicsIDs[i];
    } //end for
  } //end setGInvEquationIDs

  //! Define CMS equation IDs
  /*! Define the list of physics IDs associated with the CMS equations in ascending order.
     \param numCmsEquations (In) Number of CMS equations per node.
     \param physicsIDs (In) List of physics IDs associated with the CMS equations in ascending order.
  */
  void
  setCmsEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (cmsEquations_.size()!=0)
    {
      return; // Already been here
    } //end if
    cmsEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++)
    {
      cmsEquations_[i] = physicsIDs[i];
    } //end for
  } //end setCmsEquationIDs

  //! Define primitive density equation IDs
  /*! Define the list of physics IDs associated with the primitive density equations in ascending order.
     \param numDensityEquations (In) Number of primitive density equations per node.
     \param physicsIDs (In) List of physics IDs associated with the primitive density equations in ascending order.
  */
  void
  setDensityEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (densityEquations_.size()!=0)
    {
      return; // Already been here
    } //end if
    densityEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++)
    {
      densityEquations_[i] = physicsIDs[i];
    } //end for
  } //end setDensityEquationIDs

  //! Define Poisson equation IDs
  /*! Define the list of physics IDs associated with the Poisson equations in ascending order.
    \param numPoissonEquations (In) Number of Poisson equations per node. (should be 1)
    \param physicsIDs (In) List of physics IDs associated with the Poisson equations.
  */
  void
  setPoissonEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (poissonEquations_.size()!=0)
    {
      return; //Already been here
    } //end if
    poissonEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++)
    {
      poissonEquations_[i] = physicsIDs[i];
    } //end for
  } //end setPoissonEquationIDs

  //! Assert that field dependence on primitive densities is linear; manager will not reset values between nonlinear solves.
  /*! This method can be called to assert that the field variable dependence on primitive densities does change from one linear solve to the next.
      In this case, we can avoid filling the associated matrix coefficients.  Calling this method with "true" will cause the problem manager not reset
      the matrix coefficients for this block and to ignore any values that are submitted for entry in this block.
     \param isLinear (In) Set to true if the field dependence is linear on primitive densities.
     \warning This method can be called at any time, but should be called before the initializeValues() method is called for the second solve; By default the manager assumes that the relationship is non-linear, so values will be reset to zero and must be refilled before each linear solve.
  */
  void
  setFieldOnDensityIsLinear
  (bool isLinear)
  {
    isLinear_ = isLinear;
  } //end setFieldOnDenityLinear

  //! Method that must be called once, when all row and column maps are set.
  /*! This method constructs all of the CrsGraph objects and the lhs and rhs vectors.
   \pre All "set" methods must be called: setNodalRowMap(), setNodalColMap(), setGEquationIDs(), setGInvEquationIDs(), setCmsEquationIDs() and setDensityEquationIDs()
   \post The problem structure is finalized and cannot be changed.
  */
  void
  finalizeBlockStructure
  ();
  //@}

  //@{ \name Matrix, lhs and rhs value setup methods for piece-wise construction.

  //! Method that must be called each time \e prior to starting matrix, lhs and rhs value insertion (usually once per nonlinear iteration).
  /*! This method zeros out the matrix, lhs and rhs values. */
  virtual void
  initializeProblemValues
  ();

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
  virtual void
  insertMatrixValue
  (LocalOrdinal ownedPhysicsID, LocalOrdinal ownedNode, LocalOrdinal boxPhysicsID, LocalOrdinal boxNode, Scalar value);

  //! Method that must be called each time matrix value insertion is complete (usually once per nonlinear iteration).
  virtual void
  finalizeProblemValues
  ();

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

  //! Check for inconsistencies in operators.
  /* \param verbose (In) Print the residual of inv(A11)*A11*x_random and inv(A22)*A22*x_random.

  //! Throws an exception if the residual error is "large".
  */
  void
  Check
  (bool verbose) const;

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
  //@}

  //@{ \name Miscellaneous support methods (used by the application, or by LOCA, or both)

  //! Apply global linear operator for all physics types, b = Ax
  /*! Computes b = Ax where A is the global linear operator and each x[i] and b[i] is associated with the ith physics block.
    \param x (In) Array of pointers such that x[i] is an array of Scalars of length numBoxNodes.
    \param b (Out) Array of pointers such that b[i] is an array of Scalars of length numOwnedNodes.

    \warning Note that x[i] is of length numBoxNodes and b[i] is of length numOwnedNodes.
  */
  virtual ArrayRCP<ArrayRCP<Scalar> >
  applyMatrix(const ArrayView<const ArrayView<const Scalar> >& x) const;

  //@}

protected:
  void insertRowA12();
  void insertRowA21();
  Array<LocalOrdinal> gEquations_;
  Array<LocalOrdinal> gInvEquations_;
  Array<LocalOrdinal> cmsEquations_;
  Array<LocalOrdinal> densityEquations_;
  Array<LocalOrdinal> poissonEquations_;
  Array<LocalOrdinal> isCmsEquation_;
  Array<LocalOrdinal> isDensityEquation_;
  Array<LocalOrdinal> isPoissonEquation_;
  RCP<P11TO> A11_;
  RCP<P22TO> A22_;
  RCP<MAT_P> A12_;
  RCP<MAT_P> A21_;
  RCP<GRAPH> A12Graph_;
  RCP<GRAPH> A21Graph_;
  RCP<MAT_P> A12Static_;
  RCP<MAT_P> A21Static_;
  RCP<INVOP> A22precond_;
  RCP<const MAP> block1RowMap_;
  RCP<const MAP> block2RowMap_;
  RCP<const MAP> cmsRowMap_;
  RCP<const MAP> densityRowMap_;
  RCP<const MAP> extraRowMap_;
  RCP<const MAP> poissonRowMap_;
  bool hasPoisson_;
  Array<GlobalOrdinal> physicsIdToSchurBlockId_;
  RCP<ScTO> schurOperator_;
  RCP<VEC> rhs1_;
  RCP<VEC> rhs2_;
  RCP<VEC> rhsSchur_;
  RCP<VEC> lhs1_;
  RCP<VEC> lhs2_;
  bool isLinear_;
  bool debug_;
  GlobalOrdinal curRowA12_;
  Array<GlobalOrdinal> indicesA12_;
  GlobalOrdinal curRowA21_;
  Array<GlobalOrdinal> indicesA21_;
  std::map<GlobalOrdinal, precScalar> curRowValuesA12_;
  Array<precScalar> valuesA12_;
  std::map<GlobalOrdinal, precScalar> curRowValuesA21_;
  Array<precScalar> valuesA21_;

  using BLPM::isBlockStructureSet_;
  using BLPM::numGlobalNodes_;
  using BLPM::numGlobalBoxNodes_;
  using BLPM::physicsOrdering_;
  using BLPM::numUnknownsPerNode_;
  using BLPM::checkPhysicsOrdering;
  using BLPM::solverOrdering_;
  using BLPM::numOwnedNodes_;
  using BLPM::ownedMap_;
  using BLPM::globalRowMap_;
  using BLPM::comm_;
  using BLPM::parameterList_;
  using BLPM::tpetraParameterList_;
  using BLPM::globalMatrix_;
  using BLPM::globalRhs_;
  using BLPM::globalLhs_;
  using BLPM::ownedToBoxImporter_;
  using BLPM::boxMap_;
  using BLPM::isGraphStructureSet_;
  using BLPM::isLinearProblemSet_;
  using BLPM::firstTime_;
  using BLPM::ownedToSolverGID;
  using BLPM::writeLhs;
  using BLPM::writePermutation;
  using BLPM::setLhs;
  using BLPM::getRhs;
  using BLPM::solver_;
#ifdef SUPPORTS_STRATIMIKOS
  using BLPM::lows_;
  using BLPM::thyraOp_;
  using BLPM::thyraRhs_;
  using BLPM::thyraLhs_;
  using BLPM::lowsFactory_;
#else
  using BLPM::problem_;
#endif
};
#endif /* DFT_POLYLINPROBMGR_HPP */
