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

#ifndef DFT_HARDSPHERELINPROBMGR_HPP
#define DFT_HARDSPHERELINPROBMGR_HPP

#include "Tpetra_Headers.hpp"
#include "dft_HardSphereA11_Tpetra_Operator.hpp"
#include "dft_HardSphereA22_Tpetra_Operator.hpp"
#include "dft_A22Matrix_Tpetra_Operator.hpp"
#include "dft_Schur_Tpetra_Operator.hpp"
#include "dft_TpetraBasicLinProbMgr.hpp"

//====================================================

//! dft_HardSphereLinProbMgr:  Solver manager class for Tramonto using Trilinos.
/*! The dft_HardSphereLinProbMgr class supports solver capabilities for Tramonto.

*/
template<class Scalar,class LocalOrdinal=int,class GlobalOrdinal=LocalOrdinal,
	 class Node=Kokkos::DefaultNode::DefaultNodeType>
class dft_HardSphereLinProbMgr:
  public virtual dft_BasicLinProbMgr<Scalar, LocalOrdinal, GlobalOrdinal, Node>
{
public:
  TYPEDEF(Scalar, LocalOrdinal, GlobalOrdinal, Node);

  typedef dft_BasicLinProbMgr<Scalar,LocalOrdinal,GlobalOrdinal,Node> BLPM;
#if MIXED_PREC == 1
  typedef halfScalar opScalar;
#elif MIXED_PREC == 0
  typedef Scalar opScalar;
#endif
  typedef dft_HardSphereA11_Tpetra_Operator<opScalar,LocalOrdinal,GlobalOrdinal,Node> HS11TO;
  typedef dft_HardSphereA22_Tpetra_Operator<opScalar,LocalOrdinal,GlobalOrdinal,Node> HS22TO;
  typedef dft_A22Matrix_Tpetra_Operator<opScalar,LocalOrdinal,GlobalOrdinal,Node> A22MTO;
  typedef dft_Schur_Tpetra_Operator<opScalar,LocalOrdinal,GlobalOrdinal,Node> ScTO;

  //@{ \name Constructors/destructor.
  //! dft_HardSphereLinProbMgr Constructor.
  /*! Initialize a linear problem manager for Tramonto
     \param numUnknownsPerNode (In) The number of unknowns tracked per node of the mesh.
     \param solverOptions (In) An array of ints defined in dft_solver_defs.h containing information to
			       guide and report solver status.
     \param solverParams (In) An array of doubles defined in dft_solver_defs.h containing information to
			       guide and report solver status.
     \param comm (In) MPI communicator that should be used by the solver.

     \param formSchurMatrix (In) Forces explicit construction of S = A22 - A21*inv(A11)*A12, if true.  Otherwise S is only applied implicitly.

     \param debug (In) Turns debug mode on if set to true, false by default.
  */
  /* dft_HardSphereLinProbMgr(int numUnknownsPerNode, int * solverOptions, double * solverParams, MPI_Comm comm, bool formSchurMatrix = false, bool debug = false);*/

  //! dft_HardSphereLinProbMgr Constructor.
  /*! Initialize a linear problem manager for Tramonto
     \param numUnknownsPerNode (In) The number of unknowns tracked per node of the mesh.
     \param parameterList (In) A Teuchos::ParameterList containing information to guide and report solver status.

     \param comm (In) MPI communicator that should be used by the solver.

     \param formSchurMatrix (In) Forces explicit construction of S = A22 - A21*inv(A11)*A12, if true.  Otherwise S is only applied implicitly.

     \param debug (In) Turns debug mode on if set to true, false by default.
  */
  dft_HardSphereLinProbMgr(LocalOrdinal numUnknownsPerNode, RCP<ParameterList> parameterList, RCP<const COMM> comm, bool formSchurMatrix = false, bool debug = false);

  //! dft_HardSphereLinProbMgr Destructor.
  /*! Completely deletes a dft_HardSphereLinProbMgr object.
  */
  virtual ~dft_HardSphereLinProbMgr();
  //@}

  //@{ \name Block structure setup methods

  //! Define independent non-local equation IDs
  /*! Define the list of physics IDs associated with the non-local density equations that have no dependence on other non-local densities in ascending order.
     \param numIndpNonLocalEquations (In) Number of independent non-local equations per node.
     \param physicsIDs (In) List of physics IDs associated with the independent non-local density equations in ascending order.
  */
  void
  setIndNonLocalEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (indNonLocalEquations_.size()!=0) {
      return; // Already been here
    }
    indNonLocalEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++) {
      indNonLocalEquations_[i] = physicsIDs[i];
    }
  }

  //! Define dependent non-local equation IDs
  /*! Define the list of physics IDs associated with the non-local density equations that are dependent on other non-local densities in ascending order.
     \param numDeppNonLocalEquations (In) Number of dependent non-local equations per node.
     \param physicsIDs (In) List of physics IDs associated with the dependent non-local density equations in ascending order.
  */
  void
  setDepNonLocalEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (depNonLocalEquations_.size()!=0) {
      return; // Already been here
    }
    depNonLocalEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++) {
      depNonLocalEquations_[i] = physicsIDs[i];
    }
  }

  //! Define primitive density equation IDs
  /*! Define the list of physics IDs associated with the primitive density equations in ascending order.
     \param numDensityEquations (In) Number of primitive density equations per node.
     \param physicsIDs (In) List of physics IDs associated with the primitive density equations in ascending order.
  */
  void
  setDensityEquationIDs
  (const ArrayView<const LocalOrdinal>& physicsIDs)
  {
    if (densityEquations_.size()!=0) {
      return; // Already been here
    }
    densityEquations_.resize(physicsIDs.size());
    for (LocalOrdinal i=0; i<physicsIDs.size(); i++) {
      densityEquations_[i] = physicsIDs[i];
    }
  }

  //! Method that allows assertion that A22 block will be diagonal matrix.
  /*! Sets a bool that allows assertion that the A22 block will be strictly diagonal.  By default, we assume a general A22 block.
    \param isDiagonal (In) Set to true if you know that A22 is diagonal. Otherwise, set it false, or do not call this method since the default is false.
  */
  void
  setA22BlockIsDiagonal
  (bool isA22Diagonal) {
    if (isBlockStructureSet_) {
      return; // Block structure is already set, calling this method is not effective.
    }
    isA22Diagonal_ = isA22Diagonal;
  }

  //! Method that must be called once, when all row and column maps are set.
  /*! This method constructs all of the Epetra_CrsGraph objects and the lhs and rhs vectors.
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

     \return Returns 0 if residual is "small", otherwise it returns -1.
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

    \return Returns 0 if no error, -1 of any problems with the file system.
  */
  virtual void
  writeMatrix
  (const char * filename, const char * matrixName, const char * matrixDescription) const;
  //@}

  //@{ \name Miscellaneous support methods (used by the application, or by LOCA, or both)

  //! Apply global linear operator for all physics types, b = Ax
  /*! Computes b = Ax where A is the global linear operator and each x[i] and b[i] is associated with the ith physics block.
    \param x (In) Array of pointers such that x[i] is an array of doubles of length numBoxNodes.
    \param b (Out) Array of pointers such that b[i] is an array of doubles of length numOwnedNodes.

    \warning Note that x[i] is of length numBoxNodes and b[i] is of length numOwnedNodes.
  */
  virtual ArrayRCP<ArrayRCP<Scalar> >
  applyMatrix
  (const ArrayView<const ArrayView<const Scalar> >& x) const;

  //@}

protected:

  void insertRowA12();
  void insertRowA21();
  Array<LocalOrdinal> indNonLocalEquations_;
  Array<LocalOrdinal> depNonLocalEquations_;
  Array<LocalOrdinal> densityEquations_;
  RCP<HS11TO> A11_;
  RCP<A22MTO> A22Matrix_;
  RCP<HS22TO> A22Diagonal_;
#if MIXED_PREC == 1
  RCP<MAT_H> A12_;
  RCP<MAT_H> A21_;
  RCP<INVOP_H> A22MatrixPrecond_;
  RCP<INVOP_H> A22DiagonalPrecond_;
  RCP<HAPINV> A22MatrixPrecondMixed_;
  RCP<HAPINV> A22DiagonalPrecondMixed_;
#elif MIXED_PREC == 0
  RCP<MAT> A12_;
  RCP<MAT> A21_;
  RCP<INVOP> A22MatrixPrecond_;
  RCP<INVOP> A22DiagonalPrecond_;
#endif
  RCP<const MAP> block1RowMap_;
  RCP<const MAP> block2RowMap_;
  RCP<const MAP> indNonLocalRowMap_;
  RCP<const MAP> depNonLocalRowMap_;
  Array<GlobalOrdinal> physicsIdToSchurBlockId_;
#if MIXED_PREC == 1
  RCP<ScTO> schurOperator_;
  RCP<HAPINV> schurOperatorMixed_;
  RCP<HOP> schurComplement_;
  RCP<VEC_H> rhs1Half_;
  RCP<VEC_H> rhs2Half_;
  RCP<VEC_H> rhsSchurHalf_;
  RCP<VEC_H> lhs1Half_;
  RCP<VEC_H> lhs2Half_;
#elif MIXED_PREC == 0
  RCP<ScTO> schurOperator_;
#endif
  RCP<VEC> rhs1_;
  RCP<VEC> rhs2_;
  RCP<VEC> rhsSchur_;
  RCP<VEC> lhs1_;
  RCP<VEC> lhs2_;
  bool isA22Diagonal_;
  bool formSchurMatrix_;
  bool debug_;
  GlobalOrdinal curRowA12_;
  Array<GlobalOrdinal> indicesA12_;
  int curRowA21_;
  Array<GlobalOrdinal> indicesA21_;
#if MIXED_PREC == 1
  typedef typename std::map<GlobalOrdinal, halfScalar>::iterator ITER;
  std::map<GlobalOrdinal, halfScalar> curRowValuesA12_;
  Array<halfScalar> valuesA12_;
  std::map<GlobalOrdinal, halfScalar> curRowValuesA21_;
  Array<halfScalar> valuesA21_;
#elif MIXED_PREC == 0
  typedef typename std::map<GlobalOrdinal, Scalar>::iterator ITER;
  std::map<GlobalOrdinal, Scalar> curRowValuesA12_;
  Array<Scalar> valuesA12_;
  std::map<GlobalOrdinal, Scalar> curRowValuesA21_;
  Array<Scalar> valuesA21_;
#endif

  using BLPM::parameterList_;
  using BLPM::isBlockStructureSet_;
  using BLPM::numGlobalNodes_;
  using BLPM::numGlobalBoxNodes_;
  using BLPM::numGlobalCoarsenedNodes_;
  using BLPM::numCoarsenedNodes_;
  using BLPM::numOwnedNodes_;
  using BLPM::physicsOrdering_;
  using BLPM::numUnknownsPerNode_;
  using BLPM::checkPhysicsOrdering;
  using BLPM::solverOrdering_;
  using BLPM::ownedToBoxImporter_;
  using BLPM::boxMap_;
  using BLPM::ownedMap_;
  using BLPM::isGraphStructureSet_;
  using BLPM::isLinearProblemSet_;
  using BLPM::firstTime_;
  using BLPM::globalRowMap_;
  using BLPM::globalMatrix_;
  using BLPM::globalRhs_;
  using BLPM::globalLhs_;
  using BLPM::comm_;
  using BLPM::getRhs;
  using BLPM::ownedNodeIsCoarsened_;
  using BLPM::boxNodeIsCoarsened_;
  using BLPM::coarsenedNodesMap_;
  using BLPM::solver_;
  using BLPM::problem_;
};

#endif /* DFT_HARDSPHERELINPROBMGR_HPP */
