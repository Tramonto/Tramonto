// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

#ifndef TPETRA_CRSMATRIXMIXED_DECL_HPP
#define TPETRA_CRSMATRIXMIXED_DECL_HPP

// TODO: row-wise insertion of entries in globalAssemble() may be more efficient

// TODO: add typeglobs: CrsMatrix<Scalar,typeglob>
// TODO: add template (template) parameter for nonlocal container (this will be part of typeglob)

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp_decl.hpp"

namespace Tpetra {

  //! \brief Sparse matrix that presents a compressed sparse row interface.
  /*!
   \tparam Scalar The type of the numerical entries of the matrix.
     (You can use real-valued or complex-valued types here, unlike in
     Epetra, where the scalar type is always \c double.)

   \tparam LocalOrdinal The type of local indices.  Same as the \c
     LocalOrdinal template parameter of \c Map objects used by this
     matrix.  (In Epetra, this is just \c int.)

   \tparam GlobalOrdinal The type of global indices.  Same as the \c
     GlobalOrdinal template parameter of \c Map objects used by this
     matrix.  (In Epetra, this is just \c int.  One advantage of
     Tpetra over Epetra is that you can use a 64-bit integer type here
     if you want to solve big problems.)

   \tparam Node A class implementing on-node shared-memory parallel
     operations.  It must implement the
     \ref kokkos_node_api "Kokkos Node API."
     The default \c Node type depends on your Trilinos build options.

   \tparam LocalMatOps A local sparse matrix operations class.  It
     must implement the \ref kokkos_crs_ops "Kokkos CRS Ops API."

   This class implements a distributed-memory parallel sparse matrix,
   and provides sparse matrix-vector multiply (including transpose)
   and sparse triangular solve operations.  It provides access by rows
   to the elements of the matrix, as if the local data were stored in
   compressed sparse row format.  (Implementations are _not_ required
   to store the data in this way internally.)  This class has an
   interface like that of \c Epetra_CrsMatrix, but also allows
   insertion of data into nonowned rows, much like \c Epetra_FECrsMatrix.

   <b>Local vs. Global</b>

   The distinction between local and global indices might confuse new
   Tpetra users.  _Global_ indices represent the rows and columns
   uniquely over the entire matrix, which may be distributed over
   multiple processes.  _Local_ indices are local to the process that
   owns them.  If global index G is owned by process P, then there is
   a unique local index L on process P corresponding to G.  If the
   local index L is valid on process P, then there is a unique global
   index G owned by P corresponding to the pair (L, P).  However,
   multiple processes might own the same global index (an "overlapping
   Map"), so a global index G might correspond to multiple (L, P)
   pairs.  In summary, local indices on a process correspond to matrix
   rows or columns owned by that process.

   We summarize the different between local and global indices because
   many of CrsMatrix's methods for adding, modifying, or accessing
   entries come in versions that take either local or global indices.
   The matrix itself may store indices either as local or global.  You
   should only use the method version corresponding to the current
   state of the matrix.  For example, \c getGlobalRowView() returns a
   view to the indices represented as global; it is incorrect to call
   this method if the matrix is storing indices as local.  Call the \c
   isGloballyIndexed() or \c isLocallyIndexed() methods to find out
   whether the matrix currently stores indices as local or global.

   Method methods that work with global indices only allow operations
   on indices owned by the calling process.  For example, methods that
   take a global row index expect that row to be owned by the calling
   process.  Access to nonlocal (i.e., not owned by the calling
   process) rows requires performing an explicit communication via the
   import/export capabilities of the \c CrsMatrix object; see \c
   DistObject.  However, the method \c insertGlobalValues() is an
   exception to this rule.  It allows you to add data to nonlocal
   rows.  These data are stored locally and communicated to the
   appropriate node on the next call to \c globalAssemble() or \c
   fillComplete() (the latter calls the former).
  */
  template <class Scalar,
	    class LocalOrdinal  = int,
	    class GlobalOrdinal = LocalOrdinal,
	    class Node          = Kokkos::DefaultNode::DefaultNodeType,
	    class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps >
  class CrsMatrixMixed : public virtual CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> {
  public:
    typedef Scalar                                                    scalar_type;
    typedef LocalOrdinal                                              local_ordinal_type;
    typedef GlobalOrdinal                                             global_ordinal_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::doublePrecision   double_scalar_type;
    typedef Node                                                      node_type;
    typedef Map<LocalOrdinal,GlobalOrdinal,Node>                      map_type;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> matrix_type;
    // backwards compatibility defines both of these
    typedef LocalMatOps                                              mat_vec_type;
    typedef LocalMatOps                                              mat_solve_type;

    using matrix_type::isFillComplete;

    //! @name Constructor/Destructor Methods
    //@{

    /// \brief Constructor specifying fixed number of entries for each row.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of matrix
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrixMixed (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
		    size_t maxNumEntriesPerRow,
		    ProfileType pftype = DynamicProfile,
		    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying (possibly different) number of entries in each row.
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.

    CrsMatrixMixed (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
		    const ArrayRCP<const LocalOrdinal>& NumEntriesPerRowToAlloc,
		    ProfileType pftype = DynamicProfile,
		    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and fixed number of entries for each row.
    ///
    /// The column Map will be used to filter any matrix entries
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param maxNumEntriesPerRow [in] Maximum number of matrix
    ///   entries per row.  If pftype==DynamicProfile, this is only a
    ///   hint, and you can set this to zero without affecting
    ///   correctness.  If pftype==StaticProfile, this sets the amount
    ///   of storage allocated, and you cannot exceed this number of
    ///   entries in any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrixMixed (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
		    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
		    size_t maxNumEntriesPerRow,
		    ProfileType pftype = DynamicProfile,
		    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying column Map and number of entries in each row.
    ///
    /// The column Map will be used to filter any matrix indices
    /// inserted using insertLocalValues() or insertGlobalValues().
    ///
    /// \param rowMap [in] Distribution of rows of the matrix.
    ///
    /// \param colMap [in] Distribution of columns of the matrix.
    ///
    /// \param NumEntriesPerRowToAlloc [in] Maximum number of matrix
    ///   entries to allocate for each row.  If
    ///   pftype==DynamicProfile, this is only a hint.  If
    ///   pftype==StaticProfile, this sets the amount of storage
    ///   allocated, and you cannot exceed the allocated number of
    ///   entries for any row.
    ///
    /// \param pftype [in] Whether to allocate storage dynamically
    ///   (DynamicProfile) or statically (StaticProfile).
    ///
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    CrsMatrixMixed (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& rowMap,
		    const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& colMap,
		    const ArrayRCP<const LocalOrdinal>& NumEntriesPerRowToAlloc,
		    ProfileType pftype = DynamicProfile,
		    const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    /// \brief Constructor specifying a previously constructed graph.
    ///
    /// Calling this constructor fixes the graph structure of the
    /// sparse matrix.  We say in this case that the matrix has a
    /// "static graph."  If you create a CrsMatrixMixed with this
    /// constructor, you are not allowed to insert new entries into
    /// the matrix, but you are allowed to change values in the
    /// matrix.
    ///
    /// The given graph must be fill complete.  Note that calling
    /// resumeFill() on the graph makes it not fill complete, even if
    /// you had previously called fillComplete() on the graph.  In
    /// that case, you must call fillComplete() on the graph again
    /// before invoking this CrsMatrixMixed constructor.
    ///
    /// \param graph [in] The graph structure of the sparse matrix.
    ///   The graph <i>must</i> be fill complete.
    /// \param params [in/out] Optional list of parameters.  If not
    ///   null, any missing parameters will be filled in with their
    ///   default values.
    explicit CrsMatrixMixed (const Teuchos::RCP<const CrsGraph<LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >& graph,
			     const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Destructor.
    virtual ~CrsMatrixMixed();

    //@}
    //! @name Methods implementing Operator
    //@{

    //! \brief Computes the sparse matrix-multivector multiplication.
    /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
      - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
    */
    void applyDouble(const MultiVector<double_scalar_type,LocalOrdinal,GlobalOrdinal,Node> & X,
		     MultiVector<double_scalar_type,LocalOrdinal,GlobalOrdinal,Node> &Y,
		     Teuchos::ETransp mode = Teuchos::NO_TRANS,
		     double_scalar_type alpha = ScalarTraits<double_scalar_type>::one(),
		     double_scalar_type beta = ScalarTraits<double_scalar_type>::zero()) const;

    //@}

  protected:

    // a wrapper around multiply, for use in apply; it contains a non-owning RCP to *this, therefore, it is not allowed
    // to persist past the destruction of *this. therefore, WE MAY NOT SHARE THIS POINTER.
    RCP< const CrsMatrixMultiplyOp<double_scalar_type,Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > doubleScalarMultiplyOp_;



  }; // class CrsMatrixMixed

} // namespace Tpetra

#endif
