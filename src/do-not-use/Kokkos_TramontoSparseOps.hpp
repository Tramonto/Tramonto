//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef KOKKOS_TRAMONTOSPARSEOPS_HPP
#define KOKKOS_TRAMONTOSPARSEOPS_HPP

#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Describable.hpp>
#include <iterator>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_TramontoSparseSolveKernelOps.hpp"
#include "Kokkos_TramontoSparseMultiplyKernelOps.hpp"

namespace KokkosClassic {

  /// \class TramontoHostLocalSparseOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  /// \tparam Allocator The allocator to use when allocating sparse
  ///   matrix data.  Depending on the particular Allocator, this may
  ///   or may not do first-touch initialization.  We use first-touch
  ///   allocation by default; if you don't want this, you should use
  ///   details::DefaultCRSAllocator as the Allocator type.
  template <class Scalar,
            class Ordinal,
            class Node,
            class Allocator = details::FirstTouchCRSAllocator>
  class TramontoHostLocalSparseOps : public Teuchos::Describable {
  public:
    //! \name Typedefs and structs
    //@{

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator> sparse_ops_type;

    //! Typedef for local graph class
    template <class O, class N>
    struct graph {
      typedef DefaultCrsGraph<O,N> graph_type;
    };

    //! Typedef for local matrix class
    template <class S, class O, class N>
    struct matrix {
      typedef DefaultCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Local sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for local
    /// sparse operations for a scalar type S2, which may be different
    /// from \c Scalar.
    ///
    /// This class' typedef is used by Tpetra::CrsMatrix to bind a
    /// potentially "void" scalar type to the appropriate scalar.  The
    /// other_type typedef tells Tpetra::CrsMatrix which local sparse
    /// ops type to use, as a function of Tpetra's Scalar template
    /// parameter.
    ///
    /// Other local sparse ops implementations (especially those that
    /// wrap third-party libraries implementing sparse kernels) might
    /// use this to provide a "fall-back" sparse ops implementation of
    /// a possibly different type, if the third-party library does not
    /// support scalar type S2.
    ///
    /// In the case of TramontoHostLocalSparseOps, the other_type typedef
    /// always specifies a specialization of \c TramontoHostLocalSparseOps,
    /// regardless of the scalar type S2.  This is not necessarily
    /// true of other implementations of local sparse ops, so Tpetra
    /// developers should always get their local sparse ops type from
    /// the other_type typedef.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef TramontoHostLocalSparseOps<S2,Ordinal,Node,Allocator> other_type;
    };

    /// \brief Local sparse operations type for a different ordinal type.
    ///
    /// The bind_scalar struct defines the type responsible for local
    /// sparse operations for an ordinal type O2, which may be
    /// different from \c Ordinal.
    ///
    /// This is used by Tpetra::CrsMatrix to "bind" the local sparse
    /// ops type, given its own (Local)Ordinal type.  In the case of
    /// TramontoHostLocalSparseOps, the other_type typedef always specifies
    /// a specialization of TramontoHostLocalSparseOps, regardless of the
    /// ordinal type S2.  This is not necessarily true of other
    /// implementations of local sparse ops, so Tpetra developers
    /// should always get their local sparse ops type from the
    /// other_type typedef.
    ///
    /// Other local sparse ops implementations (especially those that
    /// wrap third-party libraries implementing sparse kernels) might
    /// use this to provide a "fall-back" sparse ops implementation of
    /// a possibly different type, if the third-party library does not
    /// support ordinal type O2.
    ///
    /// \tparam O2 An ordinal type possibly different from \c Ordinal.
    template <class O2>
    struct bind_ordinal {
      typedef TramontoHostLocalSparseOps<Scalar,O2,Node,Allocator> other_type;
    };

    //@}
    //! \name Constructors and destructor
    //@{

    /// \brief Constructor, with default parameters.
    ///
    /// We syntactically forbid setting parameters after construction,
    /// since setting parameters after calling setGraphAndMatrix()
    /// would require reorganizing the already optimized sparse matrix
    /// storage.  If you want to set nondefault values of parameters,
    /// you must use the constructor that takes a ParameterList.
    ///
    /// \param node [in/out] Kokkos Node instance.
    TramontoHostLocalSparseOps(const RCP<Node> &node);

    /// \brief Constructor, with custom parameters.
    ///
    /// Both this constructor and finalizeGraphAndMatrix() accept a
    /// ParameterList.  However, those sets of parameters are
    /// different.  The constructor's parameters concern the
    /// algorithm, and the parameters for finalizeGraphAndMatrix()
    /// concern the data structure.  It's possible to use different
    /// algorithms with the same data structure.
    ///
    /// \param node [in/out] Kokkos Node instance.
    ///
    /// \param params [in/out] Parameters for the solve.  We fill in
    ///   the given ParameterList with its default values, but we
    ///   don't keep it around.  (This saves a bit of memory.)
    TramontoHostLocalSparseOps(const RCP<Node> &node, Teuchos::ParameterList& params);

    /// \brief "Sum constructor": compute *this = alpha*A + beta*B.
    ///
    /// The resulting matrix shares the Node instance and copies the
    /// parameters of the matrix A.
    TramontoHostLocalSparseOps (const Scalar& alpha,
                          const TramontoHostLocalSparseOps<Scalar, Ordinal, Node, Allocator>& A,
                          const Scalar& beta,
                          const TramontoHostLocalSparseOps<Scalar, Ordinal, Node, Allocator>& B);

    //! Destructor
    ~TramontoHostLocalSparseOps();

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os << "KokkosClassic::TramontoHostLocalSparseOps<"
         << "Scalar=" << TypeNameTraits<Scalar>::name()
         << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //! Write a possibly more verbose description of this instance to out.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
    {
      using Teuchos::EVerbosityLevel;
      using Teuchos::includesVerbLevel;
      using Teuchos::OSTab;
      using Teuchos::rcpFromRef;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;
      using std::endl;
      typedef Teuchos::ArrayRCP<const size_t>::size_type size_type;

      // Interpret the default verbosity level as VERB_MEDIUM.
      const EVerbosityLevel vl =
        (verbLevel == VERB_DEFAULT) ? VERB_MEDIUM : verbLevel;

      if (vl == VERB_NONE) {
        return;
      }
      else if (includesVerbLevel (vl, VERB_LOW)) { // vl >= VERB_LOW
        out << this->description();

        if (includesVerbLevel (vl, VERB_MEDIUM)) { // vl >= VERB_MEDIUM
          out << " {" << endl;
          OSTab tab1 (rcpFromRef (out));

          out << "isInitialized_: " << isInitialized_ << endl;
          if (isInitialized_) {
            std::string triUplo ("INVALID");
            if (tri_uplo_ == Teuchos::UNDEF_TRI) {
              triUplo = "UNDEF_TRI";
            }
            else if (tri_uplo_ == Teuchos::LOWER_TRI) {
              triUplo = "LOWER_TRI";
            }
            else if (tri_uplo_ == Teuchos::UPPER_TRI) {
              triUplo = "UPPER_TRI";
            }
            std::string unitDiag ("INVALID");
            if (unit_diag_ == Teuchos::NON_UNIT_DIAG) {
              unitDiag = "NON_UNIT_DIAG";
            }
            else if (unit_diag_ == Teuchos::UNIT_DIAG) {
              unitDiag = "UNIT_DIAG";
            }

            out << "numRows_: " << numRows_ << endl
                << "numCols_: " << numCols_ << endl
                << "isEmpty_: " << isEmpty_ << endl
                << "tri_uplo_: " << triUplo << endl
                << "unit_diag_: " << unitDiag << endl;
            if (big_ptrs_.size() > 0) {
              out << "numEntries: " << big_ptrs_[big_ptrs_.size()-1] << endl;
            }
            else if (sml_ptrs_.size() > 0) {
              out << "numEntries: " << sml_ptrs_[sml_ptrs_.size()-1] << endl;
            }
            else {
              out << "numEntries: 0" << endl;
            }

            if (includesVerbLevel (vl, VERB_EXTREME)) { // vl >= VERB_EXTREME
              // Only print out all the sparse matrix's data in
              // extreme verbosity mode.
              out << "ptrs: [";
              if (big_ptrs_.size() > 0) {
                for (size_type i = 0; i < big_ptrs_.size (); ++i) {
                  out << big_ptrs_[i];
                  if (i + 1 < big_ptrs_.size ()) {
                    out << ", ";
                  }
                }
              }
              else {
                for (size_type i = 0; i < sml_ptrs_.size (); ++i) {
                  out << sml_ptrs_[i];
                  if (i + 1 < sml_ptrs_.size ()) {
                    out << ", ";
                  }
                }
              }
              out << "]" << endl << "inds_: [";
              for (size_type i = 0; i < inds_.size (); ++i) {
                out << inds_[i];
                if (i + 1 < inds_.size ()) {
                  out << ", ";
                }
              }
              out << "]" << endl << "vals_: [";
              for (size_type i = 0; i < vals_.size (); ++i) {
                out << vals_[i];
                if (i + 1 < vals_.size ()) {
                  out << ", ";
                }
              }
              out << "]" << endl;
            } // vl >= VERB_EXTREME
          } // if is initialized
        } // vl >= VERB_MEDIUM
      } // vl >= VERB_LOW
    }

    /// \brief Convert to dense matrix and return.
    ///
    /// \warning This method is for debugging only.  It uses a lot of
    ///   memory.  Users should never call this method.  Do not rely
    ///   on this method continuing to exist in future releases.
    ///
    /// \warning SerialDenseMatrix currently requires Ordinal=int
    ///   indices for BLAS and LAPACK compatibility.  We make no
    ///   attempt to check whether Ordinal -> int conversions
    ///   overflow.
    Teuchos::RCP<Teuchos::SerialDenseMatrix<int, scalar_type> >
    asDenseMatrix () const;

    //@}
    //! \name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of graph and matrix
    //@{

    /// \brief Allocate and initialize the storage for the row pointers.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param numEntriesPerRow [in] numEntriesPerRow[i] is the number
    ///   of entries in row i, for all rows of the local sparse matrix.
    static ArrayRCP<size_t> allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &numEntriesPerRow)
    {
      return Allocator::allocRowPtrs(node,numEntriesPerRow);
    }

    /// \brief Allocate storage for column indices or matrix values.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param rowPtrs [in] The array of row offsets; the 'ptr' array
    ///   in the compressed sparse row storage format.  rowPtrs.size()
    ///   is one plus the number of rows in the local sparse matrix.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs)
    {
      return Allocator::template allocStorage<T,size_t>(node,rowPtrs);
    }

    /// \brief Finalize the graph.
    ///
    /// \param uplo [in] Whether the structure of the graph is lower
    ///   triangular (Teuchos::LOWER_TRI), upper triangular
    ///   (Teuchos::UPPER_TRI), or neither (Teuchos::UNDEF_TRI).
    /// \param diag [in] Whether the graph has an implicitly stored
    ///   diagonal (Teuchos::UNIT_DIAG) or does not
    ///   (Teuchos::NON_UNIT_DIAG).  This currently only affects
    ///   sparse triangular solve.
    /// \param graph [in/out] The graph to finalize.
    /// \param params [in/out] Parameters for finalization.
    static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params);

    /// \brief Finalize the matrix of an already-finalized graph.
    ///
    /// \param graph [in] The graph, which must have already been
    ///   finalized using finalizeGraph().
    /// \param matrix [in/out] The matrix to finalize.  It must have
    ///   been created with the above graph as its graph.
    /// \params [in/out] Parameters for finalization.
    static void finalizeMatrix(const DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    /// \brief Finalize a graph and a matrix at the same time.
    ///
    /// \param uplo [in] Whether the structure of the graph and matrix
    ///   is lower triangular (Teuchos::LOWER_TRI), upper triangular
    ///   (Teuchos::UPPER_TRI), or neither (Teuchos::UNDEF_TRI).
    /// \param diag [in] Whether the matrix has an implicitly stored
    ///   unit diagonal (Teuchos::UNIT_DIAG) or does not
    ///   (Teuchos::NON_UNIT_DIAG).  This currently only affects
    ///   sparse triangular solve.
    /// \param graph [in/out] The graph to finalize.
    /// \param matrix [in/out] The matrix to finalize.  It must have
    ///   been created with the above graph as its graph.
    /// \param params [in/out] Parameters for finalization.
    ///
    /// Both the constructor and this method accept a ParameterList.
    /// However, those sets of parameters are different.  The
    /// constructor's parameters concern the algorithm, and the
    /// parameters for this method concern the data structure.  It's
    /// possible to use different algorithms with the same data
    /// structure.
    static void finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params);

    /// \brief Initialize sparse operations with a graph and matrix.
    ///
    /// This is the point at which this object initializes its
    /// internal representation of the local sparse matrix.  In some
    /// cases, this merely involves asking the graph and matrix for
    /// pointers to their data.  In other cases, this involves copying
    /// the data into a completely different storage format.  Whatever
    /// happens is an implementation detail of this object.
    void setGraphAndMatrix(const RCP<const DefaultCrsGraph<Ordinal,Node> > &graph,
                           const RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &matrix);

    //@}
    //! @name Computational methods
    //@{

    /// \brief Y := alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, overwriting Y with the result.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \note This method does not respect the implicit unit diagonal
    ///   indication.  If you want to simulate having an implicitly
    ///   stored unit diagonal for the operation Y := A*X, you must
    ///   compute Y := X + A*X instead.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [out] Result multivector. Contents will be overwritten.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Y := beta * Y + alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, accumulating the result into Y.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \note This method does not respect the implicit unit diagonal
    ///   indication.  If you want to simulate having an implicitly
    ///   stored unit diagonal for the operation Y := A*X, you must
    ///   compute Y := X + A*X instead.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param beta [in] Scalar constant \f$\beta\f$ by which to
    ///   multiply Y when summing with the result of the sparse
    ///   matrix-(multi)vector multiply.
    ///
    /// \param Y [in/out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              RangeScalar beta,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Solve Y = Op(A) X for X, where we assume A is triangular.
    ///
    /// Solve the (upper or lower) triangular system Y = Op(A) X.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to solve with the matrix, its
    ///   transpose, or its conjugate transpose (if applicable).
    ///
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    /// \brief Gauss-Seidel or SOR on \f$B = A X\f$.
    ///
    /// Apply a forward or backward sweep of Gauss-Seidel or
    /// Successive Over-Relaxation (SOR) to the linear system(s) \f$B
    /// = A X\f$.  For Gauss-Seidel, set the damping factor \c omega
    /// to 1.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in B.
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector B.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param B [in] Right-hand side(s).
    /// \param X [in/out] On input: initial guess(es).  On output:
    ///   result multivector(s).
    /// \param D [in] Inverse of diagonal entries of the matrix A.
    /// \param omega [in] SOR damping factor.  omega = 1 results in
    ///   Gauss-Seidel.
    /// \param direction [in] Sweep direction: Forward or Backward.
    ///   If you want a symmetric sweep, call this method twice, first
    ///   with direction = Forward then with direction = Backward.
    ///
    /// \note We don't include a separate "Symmetric" direction mode
    ///   in order to avoid confusion when using this method to
    ///   implement "hybrid" Jacobi + symmetric (Gauss-Seidel or SOR)
    ///   for a matrix distributed over multiple processes.  ("Hybrid"
    ///   means "Gauss-Seidel or SOR within the process, Jacobi
    ///   outside.")  In that case, interprocess communication (a
    ///   boundary exchange) must occur before both the forward sweep
    ///   and the backward sweep, so we would need to invoke the
    ///   kernel once per sweep direction anyway.
    template <class DomainScalar, class RangeScalar>
    void
    gaussSeidel (const MultiVector<DomainScalar,Node> &B,
                 MultiVector< RangeScalar,Node> &X,
                 const MultiVector<Scalar,Node> &D,
                 const RangeScalar& omega = Teuchos::ScalarTraits<RangeScalar>::one(),
                 const enum ESweepDirection direction = Forward) const;

    /// \brief "Add in place": compute <tt>*this = alpha*A + beta*(*this)</tt>.
    ///
    /// This method may choose to reuse storage of <tt>*this</tt>.
    void
    addInPlace (const Scalar& alpha,
                const TramontoHostLocalSparseOps<Scalar, Ordinal, Node, Allocator>& A,
                const Scalar& beta);
    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    TramontoHostLocalSparseOps(const TramontoHostLocalSparseOps& source);

    void getOffsets(ArrayRCP<const size_t> &ptrs) const {
      ptrs = big_ptrs_;
    }

    void getOffsets(ArrayRCP<const Ordinal> &ptrs) const {
      ptrs = sml_ptrs_;
    }

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void solvePrivate(Teuchos::ETransp trans,
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector< RangeScalar,Node> &X) const;

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void gaussSeidelPrivate (MultiVector< RangeScalar,Node> &X,
                             const MultiVector<DomainScalar,Node> &B,
                             const MultiVector<Scalar,Node> &D,
                             const RangeScalar& omega = Teuchos::ScalarTraits<RangeScalar>::one(),
                             const ESweepDirection direction = Forward) const;

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void multiplyPrivate(Teuchos::ETransp trans,
                         RangeScalar alpha,
                         const MultiVector<DomainScalar,Node> &X,
                               MultiVector<RangeScalar,Node> &Y) const;

    template <class DomainScalar, class RangeScalar, class OffsetType>
    void multiplyPrivate(Teuchos::ETransp trans,
                         RangeScalar alpha,
                         const MultiVector<DomainScalar,Node> &X,
                         RangeScalar beta,
                         MultiVector<RangeScalar,Node> &Y) const;

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    // packed CRS: array of row pointers, array of indices, array of values
    // pointers are EITHER size_t or Ordinal; one will be null
    ArrayRCP<const Ordinal> inds_;
    ArrayRCP<const size_t>  big_ptrs_;
    ArrayRCP<const Ordinal> sml_ptrs_;
    ArrayRCP<const Scalar>  vals_;

    Teuchos::EUplo  tri_uplo_;
    Teuchos::EDiag unit_diag_;

    Ordinal numRows_;
    Ordinal numCols_;
    bool isInitialized_;
    bool isEmpty_;
  };


  template <class Scalar, class Ordinal, class Node, class Allocator>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params)
  {
    graph.setMatDesc(uplo,diag);
    std::string FuncName("KokkosClassic::TramontoHostLocalSparseOps::finalizeGraph(graph,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
    // determine how many non-zeros, so that we can decide whether to reduce the offset pointer type
    ArrayRCP<const size_t> bigptrs = graph.getPointers();
    const size_t numrows = bigptrs.size() - 1,
                   numnz = bigptrs[numrows];
    if (numnz < (size_t)Teuchos::OrdinalTraits<Ordinal>::max()) {
      ArrayRCP<Ordinal> smallptrs = arcp<Ordinal>(numrows+1);
      std::copy( bigptrs.begin(), bigptrs.end(), smallptrs.begin() );
      graph.setSmallPointers(smallptrs);
    }
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeMatrix(const DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params)
  {
    // nothing much to do here
    std::string FuncName("KokkosClassic::TramontoHostLocalSparseOps::finalizeMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false,
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::finalizeGraphAndMatrix(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, DefaultCrsMatrix<Scalar,Ordinal,Node> &matrix, const RCP<ParameterList> &params)
  {
    // finalize them individually; no benefit to doing them together
    finalizeGraph(uplo,diag,graph,params);
    finalizeMatrix(graph,matrix,params);
  }


  template<class Scalar, class Ordinal, class Node, class Allocator>
  TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::TramontoHostLocalSparseOps(const RCP<Node> &node)
  : node_(node)
  , numRows_(0)
  , numCols_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize TramontoHostLocalSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::
  TramontoHostLocalSparseOps (const RCP<Node> &node, Teuchos::ParameterList& params)
  : node_(node)
  , numRows_(0)
  , numCols_(0)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    (void) params; // Not using this yet.

    // Make sure that users only specialize TramontoHostLocalSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::
  TramontoHostLocalSparseOps (const Scalar& alpha,
                        const TramontoHostLocalSparseOps<Scalar, Ordinal, Node, Allocator>& A,
                        const Scalar& beta,
                        const TramontoHostLocalSparseOps<Scalar, Ordinal, Node, Allocator>& B)
  {
    // Make sure that users only specialize TramontoHostLocalSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;

    (void) alpha;
    (void) A;
    (void) beta;
    (void) B;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "KokkosClassic::"
      "TramontoHostLocalSparseOps: sum constructor not implemented.");
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  void
  TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::
  addInPlace (const Scalar& alpha,
              const TramontoHostLocalSparseOps<Scalar, Ordinal, Node, Allocator>& A,
              const Scalar& beta)
  {
    (void) alpha;
    (void) A;
    (void) beta;
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "KokkosClassic::"
      "TramontoHostLocalSparseOps::addInPlace: Not implemented.");
  }

  template<class Scalar, class Ordinal, class Node, class Allocator>
  TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::~TramontoHostLocalSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  RCP<Node> TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::getNode() const {
    return node_;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::setGraphAndMatrix(
              const RCP<const DefaultCrsGraph<Ordinal,Node> > &opgraph,
              const RCP<const DefaultCrsMatrix<Scalar,Ordinal,Node> > &opmatrix)
  {
    std::string tfecfFuncName("setGraphAndMatrix(uplo,diag,graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " operators already initialized.");
    numRows_ = opgraph->getNumRows();
    numCols_ = opgraph->getNumCols();
    if (opgraph->isEmpty() || numRows_ == 0 || numCols_ == 0) {
      isEmpty_ = true;
    }
    else {
      isEmpty_ = false;
      big_ptrs_ = opgraph->getPointers();
      sml_ptrs_ = opgraph->getSmallPointers();
      inds_ = opgraph->getIndices();
      vals_ = opmatrix->getValues();
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        big_ptrs_ != null && sml_ptrs_ != null,
        std::logic_error, " Internal logic error: graph has small and big pointers. Please notify Kokkos team."
      )
      const size_t lenptrs = (big_ptrs_ != null ? big_ptrs_.size() : sml_ptrs_.size());
      // these checks just about the most that we can perform
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          lenptrs != (size_t)numRows_+1 || inds_.size() != vals_.size(),
          std::runtime_error, " matrix and graph seem incongruent.");
    }
    opgraph->getMatDesc( tri_uplo_, unit_diag_ );
    isInitialized_ = true;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::solvePrivate(Teuchos::ETransp trans,
                                                        const MultiVector<DomainScalar,Node> &Y,
                                                              MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    typedef TramontoSparseSolveOp<         Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar>   Op;
    typedef TramontoSparseTransposeSolveOp<Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar>  TOp;
    ArrayRCP<const OffsetType> ptrs;
    getOffsets(ptrs);
    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(unit_diag_ != Teuchos::UNIT_DIAG, std::runtime_error,
          " solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.offs    = rbh.template addConstBuffer< OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (unit_diag_ == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = ( tri_uplo_ == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        // no parallel for you
        wdp.execute();
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        wdp.offs    = rbh.template addConstBuffer< OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<    Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<     Scalar>(vals_);
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (unit_diag_ == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = ( tri_uplo_ == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::
  gaussSeidelPrivate (MultiVector<RangeScalar,Node> &X,
                      const MultiVector<DomainScalar,Node> &B,
                      const MultiVector<Scalar,Node> &D,
                      const RangeScalar& omega,
                      const ESweepDirection direction) const
  {
    typedef Teuchos::ScalarTraits<RangeScalar> STS;

    if (numRows_ == 0) {
      return; // Nothing to do.
    }
    const size_t numCols = B.getNumCols ();
    if (numCols == 0) {
      return; // Nothing to do.
    }
    if (isEmpty_) {
      // All the off-diagonal entries of A are zero, and all the
      // diagonal entries are (implicitly) 1.  Therefore compute:
      //
      // X := (1 - omega) X + omega B.
      //
      // FIXME (mfh 24 Dec 2012) This only works if RangeScalar ==
      // DomainScalar.  This is because DefaultArithmetic only has one
      // template parameter, namely one multivector type.
      typedef DefaultArithmetic<MultiVector<RangeScalar,Node> > MVT;
      MVT::GESUM (X, omega, B, (STS::one() - omega));
      return;
    }

    // Get the raw pointers to all the arrays.
    ArrayRCP<const OffsetType> ptr_wrapped;
    getOffsets (ptr_wrapped);
    const OffsetType* const ptr = ptr_wrapped.getRawPtr ();
    const Ordinal* const ind = inds_.getRawPtr ();
    const Scalar* const val = vals_.getRawPtr ();
    const DomainScalar* const b = B.getValues ().getRawPtr ();
    const size_t b_stride = B.getStride ();
    RangeScalar* const x = X.getValuesNonConst ().getRawPtr ();
    const size_t x_stride = X.getStride ();
    const Scalar* const d = D.getValues ().getRawPtr ();
    const Ordinal numRows = numRows_;

    if (numCols == 1) {
      RangeScalar x_temp;

      if (direction == Forward) {
        for (Ordinal i = 0; i < numRows; ++i) {
          x_temp = Teuchos::ScalarTraits<RangeScalar>::zero ();
          for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
            const Ordinal j = ind[k];
            const Scalar A_ij = val[k];
            x_temp += A_ij * x[j];
          }
          x[i] += omega * d[i] * (b[i] - x_temp);
        }
      } else if (direction == Backward) {
        for (Ordinal i = numRows - 1; i >= 0; --i) {
          x_temp = Teuchos::ScalarTraits<RangeScalar>::zero ();
          for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
            const Ordinal j = ind[k];
            const Scalar A_ij = val[k];
            x_temp += A_ij * x[j];
          }
          x[i] += omega * d[i] * (b[i] - x_temp);
        }
      }
    }
    else { // numCols > 1
      // mfh 20 Dec 2012: If Gauss-Seidel for multivectors with
      // multiple columns becomes important, we can add unrolled
      // implementations.  The implementation below is not unrolled.
      // It may also be reasonable to parallelize over right-hand
      // sides, if there are enough of them, especially if the matrix
      // fits in cache.
      Teuchos::Array<RangeScalar> temp (numCols);
      RangeScalar* const x_temp = temp.getRawPtr ();

      if (direction == Forward) {
        for (Ordinal i = 0; i < numRows; ++i) {
          for (size_t c = 0; c < numCols; ++c) {
            x_temp[c] = Teuchos::ScalarTraits<RangeScalar>::zero ();
          }
          for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
            const Ordinal j = ind[k];
            const Scalar A_ij = val[k];
            for (size_t c = 0; c < numCols; ++c) {
              x_temp[c] += A_ij * x[j + x_stride*c];
            }
          }
          for (size_t c = 0; c < numCols; ++c) {
            x[i + x_stride*c] += omega * d[i] * (b[i + b_stride*c] - x_temp[c]);
          }
        }
      } else if (direction == Backward) { // backward mode
        for (Ordinal i = numRows - 1; i >= 0; --i) {
          for (size_t c = 0; c < numCols; ++c) {
            x_temp[c] = Teuchos::ScalarTraits<RangeScalar>::zero ();
          }
          for (OffsetType k = ptr[i]; k < ptr[i+1]; ++k) {
            const Ordinal j = ind[k];
            const Scalar A_ij = val[k];
            for (size_t c = 0; c < numCols; ++c) {
              x_temp[c] += A_ij * x[j + x_stride*c];
            }
          }
          for (size_t c = 0; c < numCols; ++c) {
            x[i + x_stride*c] += omega * d[i] * (b[i + b_stride*c] - x_temp[c]);
          }
        }
      }
    }
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::solve(Teuchos::ETransp trans,
                                                        const MultiVector<DomainScalar,Node> &Y,
                                                              MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " this solve was not fully initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumCols() != (size_t)Y.getNumCols(),
        std::runtime_error, " Left hand side and right hand side multivectors have differing numbers of vectors."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)X.getNumRows() < (size_t)numRows_,
        std::runtime_error, " Left-hand-side multivector does not have enough rows. "
                            "Likely cause is that the column map was not provided to "
                            "the Tpetra::CrsMatrix in the case of an implicit unit diagonal."
    );
    if (big_ptrs_ != null) {
      solvePrivate<DomainScalar,RangeScalar,size_t>(trans,Y,X);
    }
    else {
      solvePrivate<DomainScalar,RangeScalar,Ordinal>(trans,Y,X);
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void
  TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::
  gaussSeidel (const MultiVector<DomainScalar,Node> &B,
               MultiVector< RangeScalar,Node> &X,
               const MultiVector<Scalar,Node> &D,
               const RangeScalar& dampingFactor,
               const ESweepDirection direction) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      isInitialized_ == false,
      std::runtime_error,
      "KokkosClassic::TramontoHostLocalSparseOps::gaussSeidel: "
      "The solve was not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      isEmpty_ && unit_diag_ != Teuchos::UNIT_DIAG && numRows_ > 0,
      std::runtime_error,
      "KokkosClassic::TramontoHostLocalSparseOps::gaussSeidel: Local Gauss-Seidel with a "
      "sparse matrix with no entries, but a nonzero number of rows, is only "
      "valid if the matrix has an implicit unit diagonal.  This matrix does "
      "not.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      (size_t) X.getNumCols() != (size_t) B.getNumCols(),
      std::runtime_error,
      "KokkosClassic::TramontoHostLocalSparseOps::gaussSeidel: "
      "The multivectors B and X have different numbers of vectors.  "
      "X has " << X.getNumCols() << ", but B has " << B.getNumCols() << ".");
    TEUCHOS_TEST_FOR_EXCEPTION(
      (size_t) X.getNumRows() < (size_t) numRows_,
      std::runtime_error,
      "KokkosClassic::TramontoHostLocalSparseOps::gaussSeidel: "
      "The input/output multivector X does not have enough rows for the "
      "matrix.  X has " << X.getNumRows() << " rows, but the (local) matrix "
      "has " << numRows_ << " rows.  One possible cause is that the column map "
      "was not provided to the Tpetra::CrsMatrix in the case of a matrix with "
      "an implicitly stored unit diagonal.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      (size_t) B.getNumRows() < (size_t) numRows_,
      std::runtime_error,
      "KokkosClassic::TramontoHostLocalSparseOps::gaussSeidel: "
      "The input multivector B does not have enough rows for the "
      "matrix.  B has " << B.getNumRows() << " rows, but the (local) matrix "
      "has " << numRows_ << " rows.");

    if (big_ptrs_ != null) {
      gaussSeidelPrivate<DomainScalar,RangeScalar,size_t> (X, B, D,
                                                           dampingFactor,
                                                           direction);
    }
    else {
      gaussSeidelPrivate<DomainScalar,RangeScalar,Ordinal> (X, B, D,
                                                            dampingFactor,
                                                            direction);
    }
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::multiplyPrivate(
                                Teuchos::ETransp trans,
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X,
                                      MultiVector<RangeScalar ,Node> &Y) const
  {
    // the 1 template parameter below means that beta is not used in computations
    // and the output multivector enjoys overwrite semantics (i.e., will overwrite data/NaNs in Y)
    typedef TramontoSparseMultiplyOp<         Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 1>  Op;
    typedef TramontoSparseTransposeMultiplyOp<Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 1> TOp;
    ArrayRCP<const OffsetType> ptrs;
    getOffsets(ptrs);
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        node_->template parallel_for<Op>(0,numRows_,wdp);
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X,
                                      MultiVector<RangeScalar ,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, " X and Y do not have the same number of columns.");
    if (big_ptrs_ != null) {
      multiplyPrivate<DomainScalar,RangeScalar,size_t>(trans,alpha,X,Y);
    }
    else {
      multiplyPrivate<DomainScalar,RangeScalar,Ordinal>(trans,alpha,X,Y);
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar, class OffsetType>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::multiplyPrivate(
                                Teuchos::ETransp trans,
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    // the 0 template parameter below means that beta is used in computations
    // and the output multivector enjoys accumulation semantics (i.e., will not overwrite data/NaNs in Y)
    typedef TramontoSparseMultiplyOp<         Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 0>  Op;
    typedef TramontoSparseTransposeMultiplyOp<Scalar,OffsetType,Ordinal,DomainScalar,RangeScalar, 0> TOp;
    ArrayRCP<const OffsetType> ptrs;
    getOffsets(ptrs);
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y
      //   <= beta * Y
      // NOTE: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      //       however, X and Y may be of different size, and we need entries to determine how to mix those potential NaNs in X into Y
      //       therefore, the best we can do is scale Y to zero. Setting Y to zero would destroy NaNs in Y, which violates the semantics of the call.
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        node_->template parallel_for<Op>(0,numRows_,wdp);
      }
      else {
        TOp wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.offs    = rbh.template addConstBuffer<  OffsetType>(ptrs);
        wdp.inds    = rbh.template addConstBuffer<     Ordinal>(inds_);
        wdp.vals    = rbh.template addConstBuffer<      Scalar>(vals_);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        wdp.numRHS  = X.getNumCols();
        rbh.end();
        // no parallel for you
        wdp.execute();
      }
    }
    return;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == false,
        std::runtime_error, " sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        X.getNumCols() != Y.getNumCols(),
        std::runtime_error, " X and Y do not have the same number of columns.");
    if (big_ptrs_ != null) {
      multiplyPrivate<DomainScalar,RangeScalar,size_t>(trans,alpha,X,beta,Y);
    }
    else {
      multiplyPrivate<DomainScalar,RangeScalar,Ordinal>(trans,alpha,X,beta,Y);
    }
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  Teuchos::RCP<Teuchos::SerialDenseMatrix<int, typename TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::scalar_type> >
  TramontoHostLocalSparseOps<Scalar,Ordinal,Node,Allocator>::
  asDenseMatrix () const
  {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    typedef Teuchos::OrdinalTraits<ordinal_type> OTO;
    typedef Teuchos::ScalarTraits<scalar_type> STS;
    typedef Teuchos::SerialDenseMatrix<int, scalar_type> dense_matrix_type;

    RCP<dense_matrix_type> A_ptr =
      rcp (new dense_matrix_type (numRows_, numCols_));
    dense_matrix_type& A = *A_ptr; // for notational convenience

    if (big_ptrs_.size() > 0) {
      ArrayRCP<const size_t> ptr = big_ptrs_;
      for (ordinal_type i = OTO::zero(); i < numRows_; ++i) {
        for (size_t k = ptr[i]; k < ptr[i+1]; ++k) {
          const ordinal_type j = inds_[k];
          const scalar_type A_ij = vals_[k];
          A(i,j) += A_ij;
        }
        if (unit_diag_ == Teuchos::UNIT_DIAG) {
          // Respect whatever is in the sparse matrix, even if it is wrong.
          // This is helpful for debugging.
          A(i,i) += STS::one ();
        }
      }
    }
    else if (sml_ptrs_.size() > 0) {
      ArrayRCP<const ordinal_type> ptr = sml_ptrs_;
      for (ordinal_type i = OTO::zero(); i < numRows_; ++i) {
        for (ordinal_type k = ptr[i]; k < ptr[i+1]; ++k) {
          const ordinal_type j = inds_[k];
          const scalar_type A_ij = vals_[k];
          A(i,j) += A_ij;
        }
        if (unit_diag_ == Teuchos::UNIT_DIAG) {
          // Respect whatever is in the sparse matrix, even if it is wrong.
          // This is helpful for debugging.
          A(i,i) += STS::one ();
        }
      }
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "KokkosClassic::DefaultHost"
        "LocalSparseOps::asDenseMatrix: both big_ptrs_ and sml_ptrs_ are empty.  "
        "Please report this bug to the Kokkos developers.");
    }
    return A_ptr;
  }


  /// \brief Partial specialization of TramontoHostLocalSparseOps for Scalar=void.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  /// \tparam Allocator The allocator type.
  ///
  /// \warning This partial specialization is _not_ for users.  Kokkos
  ///   developers should see the discussion below explaining why we
  ///   need this partial specialization.
  ///
  /// Developer documentation
  /// =======================
  ///
  /// We include a partial specialization as a work-around for a
  /// Windows MSVC compilation problem reported on 08 Aug 2012 by
  /// Brent Perschbacher.  The issue is that MSVC is attempting to
  /// compile the generic methods for Scalar=void, since we do refer
  /// to the type for Scalar=void in e.g., Tpetra::CrsGraph.  However,
  /// whenever we refer to the Scalar=void case, we only reference the
  /// typedefs and inner classes inside, not the methods.  Other
  /// compilers do not attempt to compile methods of a template class
  /// that aren't called; MSVC apparently does.
  ///
  /// Kokkos developers must imitate TramontoHostLocalSparseOps by providing
  /// their own partial specializations of their local sparse kernels
  /// classes for the Scalar=void case.
  ///
  /// gcc 4.5.1 says that "default template arguments may not be used
  /// in partial specializations," so we aren't allowed to specify a
  /// default Allocator.
  template <class Ordinal, class Node, class Allocator>
  class TramontoHostLocalSparseOps<void, Ordinal, Node, Allocator> : public Teuchos::Describable {
  public:
    //! \name Typedefs and structs
    //@{

    //! The type of the individual entries of the sparse matrix.
    typedef void scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node node_type;
    //! The type of this object, the sparse operator object
    typedef TramontoHostLocalSparseOps<void, Ordinal, Node, Allocator> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef DefaultCrsGraph<O,N> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef DefaultCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// scalar_type.  (In fact, it _should_ be different than
    /// scalar_type=void in this case.  The intended use case is to
    /// rebind scalar_type=void to the different scalar type S2 of a
    /// matrix.)
    ///
    /// This always specifies a specialization of \c
    /// TramontoHostLocalSparseOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c scalar_type.
    template <class S2>
    struct bind_scalar {
      typedef TramontoHostLocalSparseOps<S2,Ordinal,Node,Allocator> other_type;
    };

    /// \brief Sparse operations type for a different ordinal type.
    ///
    /// The bind_ordinal struct defines the type responsible for
    /// sparse operations for an ordinal type O2, which may be
    /// different from Ordinal.
    ///
    /// This always specifies a specialization of \c
    /// TramontoHostLocalSparseOps, regardless of the ordinal type O2.
    ///
    /// \tparam O2 An ordinal type possibly different from \c Ordinal.
    template <class O2>
    struct bind_ordinal {
      typedef TramontoHostLocalSparseOps<void, O2, Node, Allocator> other_type;
    };

    //@}
    //! \name Constructors/Destructor
    //@{

    //! Constructor that takes a Kokkos Node: DO NOT CALL (Scalar=void specialization).
    TramontoHostLocalSparseOps (const RCP<Node> &node) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Someone attempted to "
        "instantiate KokkosClassic::TramontoHostLocalSparseOps with Scalar=void.  "
        "This is not allowed.  "
        "The Scalar=void specialization exists only for its typedefs.  "
        "Please report this bug to the Kokkos developers.");
    }

    //! Constructor that takes a Kokkos Node and parameters: DO NOT CALL (Scalar=void specialization).
    TramontoHostLocalSparseOps (const RCP<Node> &node, Teuchos::ParameterList& params) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Someone attempted to "
        "instantiate KokkosClassic::TramontoHostLocalSparseOps with Scalar=void.  "
        "This is not allowed.  "
        "The Scalar=void specialization exists only for its typedefs.  "
        "Please report this bug to the Kokkos developers.");
    }

    //! Destructor.
    ~TramontoHostLocalSparseOps() {
      // We don't throw an exception here, because throwing exceptions
      // in a destructor may cause the application to terminate [1].
      // However, it's impossible that execution will reach this
      // point, since all the constructors throw exceptions.
      //
      // [1] http://www.parashift.com/c++-faq/dtors-shouldnt-throw.html
    }

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os <<  "KokkosClassic::TramontoHostLocalSparseOps<"
         << "Scalar=void"
         << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //! Write a possibly more verbose description of this instance to out.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
    {
      using Teuchos::includesVerbLevel;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;

      // Interpret the default verbosity level as VERB_MEDIUM.
      const Teuchos::EVerbosityLevel vl =
        (verbLevel == VERB_DEFAULT) ? VERB_MEDIUM : verbLevel;

      if (vl == VERB_NONE) {
        return;
      }
      else if (includesVerbLevel (vl, VERB_LOW)) { // vl >= VERB_LOW
        out << this->description() << std::endl;
      }
    }

    //@}
    //! \name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const {
      // You're not supposed to instantiate this object, so we always
      // return null here.
      return Teuchos::null;
    }

    //@}
    //! @name Initialization of graph and matrix
    //@{

    /// \brief Allocate and initialize the storage for the row offsets.
    ///
    /// \note This is still implemented in the Scalar=void
    ///   specialization, since Tpetra::CrsGraph may use it for
    ///   allocating its row offsets.  Since it's a class method, we
    ///   may call it without needing to instantiate a
    ///   TramontoHostLocalSparseOps instance.
    static ArrayRCP<size_t> allocRowPtrs(const RCP<Node> &node, const ArrayView<const size_t> &numEntriesPerRow)
    {
      return Allocator::allocRowPtrs(node,numEntriesPerRow);
    }

    /// \brief Allocate and initialize the storage for graph or matrix storage.
    ///
    /// \note This is still implemented in the Scalar=void
    ///   specialization, since Tpetra::CrsGraph may use it for
    ///   allocating the column indices (T=Ordinal).  Since it's a
    ///   class method, we may call it without needing to instantiate
    ///   a TramontoHostLocalSparseOps instance.
    template <class T>
    static ArrayRCP<T> allocStorage(const RCP<Node> &node, const ArrayView<const size_t> &rowPtrs)
    {
      return Allocator::template allocStorage<T,size_t>(node,rowPtrs);
    }

    /// \brief Finalize a graph.
    ///
    /// \note This is still implemented in the Scalar=void
    ///   specialization, since Tpetra::CrsGraph may use it for
    ///   finalizing the graph structure.  Since it's a class method,
    ///   we may call it without needing to instantiate a
    ///   TramontoHostLocalSparseOps instance.
    static void finalizeGraph(Teuchos::EUplo uplo, Teuchos::EDiag diag, DefaultCrsGraph<Ordinal,Node> &graph, const RCP<ParameterList> &params) {
      using Teuchos::ArrayRCP;
      using Teuchos::as;

      std::string FuncName("KokkosClassic::TramontoHostLocalSparseOps::finalizeGraph(graph,params)");

      graph.setMatDesc(uplo,diag);
      TEUCHOS_TEST_FOR_EXCEPTION(graph.isInitialized() == false,
        std::runtime_error, FuncName << ": graph has not yet been initialized.");

      // determine how many non-zeros, so that we can decide whether to reduce the offset pointer type
      ArrayRCP<const size_t> bigptrs = graph.getPointers();
      const size_t numrows = bigptrs.size() - 1,
        numnz = bigptrs[numrows];
      if (numnz < as<size_t> (Teuchos::OrdinalTraits<Ordinal>::max())) {
        ArrayRCP<Ordinal> smallptrs = arcp<Ordinal>(numrows+1);
        std::copy( bigptrs.begin(), bigptrs.end(), smallptrs.begin() );
        graph.setSmallPointers(smallptrs);
      }
    }

    //@}

  private:
    //! Copy constructor (protected and unimplemented)
    TramontoHostLocalSparseOps(const TramontoHostLocalSparseOps& source);
  };

  /** \example CrsMatrix_DefaultMultiplyTests.hpp
    * This is an example that unit tests and demonstrates the implementation requirements for the DefaultLocalSparseOps class.
    */

} // namespace KokkosClassic

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */

