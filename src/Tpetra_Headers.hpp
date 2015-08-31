#ifndef TPETRA_HEADERS_H
#define TPETRA_HEADERS_H

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_MultiVector.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrixMultiplyOp.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Operator.hpp>
#include "Tpetra_OperatorApplyInverse.hpp"
#include "Tpetra_InvOperator.hpp"
#include "Tpetra_MixedOperator.hpp"
#include "Tpetra_ApplyOp.hpp"
#include "Tpetra_ParameterListConverter.hpp"
#include "Tpetra_MultiVectorConverter.hpp"
#include "Tpetra_ScalingCrsMatrix.hpp"
//#include "Kokkos_DefaultKernels.hpp"
//#include "Kokkos_TramontoSparseOps.hpp"
//#include "Kokkos_TramontoSparseMultiplyKernelOps.hpp"
//#include "Kokkos_TramontoSparseSolveKernelOps.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_VerboseObject.hpp>

#ifdef SUPPORTS_STRATIMIKOS
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_PreconditionerFactoryHelpers.hpp>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#endif

#include <BelosLinearProblem.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>

#include <Ifpack2_ILUT.hpp>
#include <Ifpack2_Diagonal.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_AdditiveSchwarz.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ParameterList;
using Teuchos::VerboseObjectBase;
using Teuchos::FancyOStream;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::rcp_dynamic_cast;

#ifdef SUPPORTS_STRATIMIKOS
using Thyra::SolveStatus;
using Thyra::createVector;
using Thyra::createMultiVector;
using Thyra::createLinearOp;
using Thyra::linearOpWithSolve;
using Thyra::prec;

using Stratimikos::DefaultLinearSolverBuilder;
#endif

//#define KDEBUG

//
// Select local ordinal and global ordinal types
//
#define L_OR                     int
#define G_OR                     int

//
// Select scalar type of matrix data and solver
//
#if LINSOLVE_PREC == 0
// Solver precision
#define SCALAR                 float
// Matrix precision
#define MAT_SCALAR             SCALAR

#elif LINSOLVE_PREC == 1
// Solver precision
#define SCALAR                 double
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR             float
#else
#define MAT_SCALAR             SCALAR
#endif

#elif LINSOLVE_PREC == 2
#include <qd/dd_real.h>
// Solver precision
#define SCALAR                 dd_real
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR             double
#else
#define MAT_SCALAR             SCALAR
#endif

#elif LINSOLVE_PREC == 3
#include <qd/qd_real.h>
// Solver precision
#define SCALAR                 qd_real
// Matrix precision
#if MIXED_PREC == 1
#define MAT_SCALAR             dd_real
#else
#define MAT_SCALAR             SCALAR
#endif
#endif

//
// Select node type
//
#if NODE_TYPE == 0
#define KOKKOS_NODE_STRING Kokkos::TPINode
#elif NODE_TYPE == 1
#define KOKKOS_NODE_STRING Kokkos::TBBNode
#elif NODE_TYPE == 2
#define KOKKOS_NODE_STRING Kokkos::OpenMP
#else
#define KOKKOS_NODE_STRING Kokkos::Serial
#endif

typedef Kokkos::Compat::KokkosDeviceWrapperNode<KOKKOS_NODE_STRING> KokkosNode;
#define NODE_STRING KokkosNode

//
// Select platform type
//
#if PLATFORM_TYPE == 0
#define PLATFORM_STRING Tpetra::SerialPlatform
#else
#define PLATFORM_STRING Tpetra::MpiPlatform
#endif

//
// Select matrix type
//
#define MATRIX_STRING Tpetra::CrsMatrix

///
// Macros for typedefs and instantiations
//
#define TRAMONTO_NODE_TYPEDEF(CLASSNAME,NAME)	\
  typedef CLASSNAME NAME;

#define TRAMONTO_PLATFORM_TYPEDEF(CLASSNAME,NODE,NAME)	\
  typedef CLASSNAME<NODE> NAME;

#define TRAMONTO_SPARSEOPS_TYPEDEF(CLASSNAME,S,LO,NODE,NAME)	\
  typedef CLASSNAME<S,LO,NODE> NAME;

#define TRAMONTO_MATRIX_TYPEDEF(CLASSNAME,S,LO,GO,NODE,NAME)	\
  typedef CLASSNAME<S,LO,GO,NODE> NAME;

#define TRAMONTO_CLASS_TYPEDEF(CLASSNAME,S,MAT,NAME)	\
  typedef CLASSNAME<S,MAT> NAME;

#define TRAMONTO_CLASS_INST(CLASSNAME,S,MAT)	\
  template class CLASSNAME<S,MAT>;

#define TRAMONTO_TYPEDEF_HELPER(CLASSNAME,NAME)	\
  TRAMONTO_NODE_TYPEDEF(NODE_STRING,NODE) \
  TRAMONTO_PLATFORM_TYPEDEF(PLATFORM_STRING,NODE,PLATFORM) \
  TRAMONTO_MATRIX_TYPEDEF(MATRIX_STRING,MAT_SCALAR,L_OR,G_OR,NODE,MAT) \
  TRAMONTO_CLASS_TYPEDEF(CLASSNAME, SCALAR, MAT, NAME)

#define TRAMONTO_INST_HELPER(CLASSNAME)	\
  TRAMONTO_NODE_TYPEDEF(NODE_STRING,NODE) \
  TRAMONTO_PLATFORM_TYPEDEF(PLATFORM_STRING,NODE,PLATFORM) \
  TRAMONTO_MATRIX_TYPEDEF(MATRIX_STRING,MAT_SCALAR,L_OR,G_OR,NODE,MAT) \
  TRAMONTO_CLASS_INST(CLASSNAME, SCALAR, MAT)

#ifdef SUPPORTS_STRATIMIKOS
#define TYPEDEF(SCALAR, LO, GO, NODE)		\
  \
  typedef Tpetra::MultiVector<SCALAR,LO,GO,Node> MV; \
  typedef Tpetra::Vector<SCALAR,LO,GO,Node> VEC; \
  typedef Tpetra::Operator<SCALAR,LO,GO,Node> OP; \
  typedef Tpetra::OperatorApplyInverse<SCALAR,LO,GO,Node> APINV; \
  typedef Tpetra::InvOperator<SCALAR,LO,GO,Node> INVOP; \
  typedef Teuchos::Comm<int> COMM; \
  typedef Tpetra::Map<LO,GO,Node> MAP; \
  typedef Tpetra::CrsMatrix<SCALAR,LO,GO,Node> MAT; \
  typedef Tpetra::Import<LO,GO,Node> IMP; \
  typedef Thyra::LinearOpBase<SCALAR> ThyraOP; \
  typedef Thyra::VectorBase<SCALAR> ThyraVEC; \
  typedef Thyra::LinearOpWithSolveBase<SCALAR> ThyraLOWS; \
  typedef Thyra::LinearOpWithSolveFactoryBase<SCALAR> ThyraLOWSFactory; \
  typedef Thyra::PreconditionerFactoryBase<SCALAR> ThyraPRECFactory; \
  typedef Thyra::PreconditionerBase<SCALAR> ThyraPREC;

#else

#define TYPEDEF(SC, MATRIX_TYPE)	\
  \
  typedef MATRIX_TYPE MAT; \
  typedef typename MAT::scalar_type MatScalar; \
  typedef typename MAT::scalar_type MSC; \
  typedef typename MAT::local_ordinal_type LocalOrdinal; \
  typedef typename MAT::local_ordinal_type LO; \
  typedef typename MAT::global_ordinal_type GlobalOrdinal; \
  typedef typename MAT::global_ordinal_type GO; \
  typedef typename MAT::node_type Node; \
  typedef typename MAT::node_type NO; \
  typedef Tpetra::RowMatrix<SC,LO,GO,NO> ROWMAT; \
  typedef Tpetra::MultiVector<SC,LO,GO,NO> MV;   \
  typedef Tpetra::Vector<SC,LO,GO,NO> VEC; \
  typedef Tpetra::Operator<SC,LO,GO,NO> OP; \
  typedef Tpetra::OperatorApplyInverse<SC,LO,GO,NO> APINV; \
  typedef Tpetra::InvOperator<SC,LO,GO,NO> INVOP; \
  typedef Tpetra::MixedOperator<SC,MSC,LO,GO,NO> MOP; \
  typedef Tpetra::CrsMatrixMultiplyOp<SC,MSC,LO,GO,NO> MMOP;        \
  typedef Teuchos::Comm<int> COMM; \
  typedef Tpetra::Map<LO,GO,NO> MAP; \
  typedef Tpetra::CrsGraph<LO,GO,NO> GRAPH;             \
  typedef Tpetra::ScalingCrsMatrix<MSC,LO,GO,NO> SCALE;     \
  typedef Tpetra::Import<LO,GO,NO> IMP; \
  typedef Belos::SolverManager<SC, MV, OP> SolMGR; \
  typedef Belos::LinearProblem<SC, MV, OP> LinPROB; \
  typedef Ifpack2::Preconditioner<SC,LO,GO,NO> PRECOND; \
  typedef Teuchos::ScalarTraits<SC> STS; \
  typedef Teuchos::ScalarTraits<MSC> STMS; \
  typedef Teuchos::OrdinalTraits<LO> OTLO; \
  typedef Teuchos::OrdinalTraits<GO> OTGO; \
  typedef typename Teuchos::ScalarTraits<SC>::halfPrecision halfScalar; \
  typedef typename Teuchos::ScalarTraits<SC>::doublePrecision doubleScalar; \
  typedef Tpetra::MultiVector<MSC,LO,GO,NO> MV_M; \
  typedef Tpetra::Vector<MSC,LO,GO,NO> VEC_M;  \
  typedef Tpetra::Vector<double,LO,GO,NO> VEC_D;	\
  typedef Tpetra::Operator<MSC,LO,GO,NO> OP_M;                          \
  typedef Ifpack2::ILUT<ROWMAT>ILUT;           \
  typedef Ifpack2::AdditiveSchwarz<ROWMAT> SCHWARZ; \
  typedef Tpetra::details::ApplyOp<SC,ILUT> ILUT_OP;                    \
  typedef Ifpack2::Diagonal<ROWMAT> DIAGONAL;   \
  typedef Tpetra::details::ApplyOp<SC,SCHWARZ> SCHWARZ_OP;	\
  typedef Tpetra::details::ApplyOp<SC,DIAGONAL> DIAGONAL_OP; \
  typedef typename std::map<GO, MSC>::iterator ITER;

#endif //SUPPORTS_STRATIMIKOS

#endif /* TPETRA_HEADERS_H */
