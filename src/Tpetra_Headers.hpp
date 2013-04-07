#ifndef TPETRA_HEADERS_H
#define TPETRA_HEADERS_H

#include "Kokkos_NodeTrace.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_MultiVector.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_OperatorApplyInverse.hpp"
#include "Tpetra_InvOperator.hpp"
#include "Tpetra_MixedOperator.hpp"
#include "Tpetra_ApplyOp.hpp"
#include "Tpetra_ParameterListConverter.hpp"
#include "Tpetra_MultiVectorConverter.hpp"
#include "Tpetra_ScalingCrsMatrix.hpp"

#include "Teuchos_Comm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_RCPDecl.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#if ENABLE_MUELU == 1
#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include "MueLu_Hierarchy.hpp"
#include <MueLu_RAPFactory.hpp>
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_Utilities.hpp"
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#include <MueLu_TpetraOperator.hpp>
#endif
#include "Kokkos_DefaultNode.hpp"
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

using Thyra::SolveStatus;
using Thyra::createVector;
using Thyra::createMultiVector;
using Thyra::createLinearOp;
using Thyra::linearOpWithSolve;
using Thyra::prec;

using Stratimikos::DefaultLinearSolverBuilder;
using Belos::ReturnType;

//#define KDEBUG
//#define SUPPORTS_STRATIMIKOS

#ifdef SUPPORTS_STRATIMIKOS
#define TYPEDEF(SCALAR, LO, GO, NODE) \
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

#define TYPEDEF(SCALAR, LO, GO, NODE) \
  \
  typedef Tpetra::MultiVector<SCALAR,LO,GO,Node> MV; \
  typedef Tpetra::Vector<SCALAR,LO,GO,Node> VEC; \
  typedef Tpetra::Operator<SCALAR,LO,GO,Node> OP; \
  typedef Tpetra::OperatorApplyInverse<SCALAR,LO,GO,Node> APINV; \
  typedef Tpetra::InvOperator<SCALAR,LO,GO,Node> INVOP; \
  typedef Tpetra::MixedOperator<SCALAR,LO,GO,Node> MOP; \
  typedef Tpetra::CrsMatrixMultiplyOp<SCALAR,SCALAR,LO,GO,Node> MMOP; \
  typedef Teuchos::Comm<int> COMM; \
  typedef Tpetra::Map<LO,GO,Node> MAP; \
  typedef Tpetra::CrsMatrix<SCALAR,LO,GO,Node> MAT; \
  typedef Tpetra::CrsGraph<LO,GO,Node> GRAPH; \
  typedef Tpetra::ScalingCrsMatrix<Scalar,LocalOrdinal,GO,Node> SCALE; \
  typedef Tpetra::Import<LO,GO,Node> IMP; \
  typedef Belos::SolverManager<SCALAR, MV, OP> SolMGR; \
  typedef Belos::LinearProblem<SCALAR, MV, OP> LinPROB; \
  typedef Ifpack2::Preconditioner<SCALAR, LO, GO, Node> PRECOND; \
  typedef typename Teuchos::ScalarTraits<SCALAR>::halfPrecision halfScalar; \
  typedef typename Teuchos::ScalarTraits<SCALAR>::doublePrecision doubleScalar; \
  typedef Tpetra::MultiVector<halfScalar,LO,GO,Node> MV_H; \
  typedef Tpetra::Vector<halfScalar,LO,GO,Node> VEC_H;		   \
  typedef Tpetra::Operator<halfScalar,LO,GO,Node> OP_H;			\
  typedef Tpetra::CrsMatrix<halfScalar,LO,GO,Node> MAT_H;			\
  typedef Tpetra::CrsMatrixMultiplyOp<SCALAR,halfScalar,LO,GO,Node> MMOP_H; \
  typedef Tpetra::ScalingCrsMatrix<halfScalar,LocalOrdinal,GO,Node> SCALE_H; \
  typedef Ifpack2::Preconditioner<halfScalar, LO, GO, Node> PRECOND_H;

#if MIXED_PREC == 1

#define TYPEDEF_MIXED(SCALAR, LO, GO, NODE) \
  \
  typedef halfScalar precScalar; \
  typedef typename std::map<GO, precScalar>::iterator ITER; \
  typedef MAT_H MAT_P; \
  typedef MMOP_H MMOP_P; \
  typedef OP_H OP_P; \
  typedef VEC_H VEC_P; \
  typedef MV_H MV_P; \
  typedef PRECOND_H PRECOND_P; \
  typedef SCALE_H SCALE_P; \
  typedef Ifpack2::ILUT<MAT_P> PRECOND_ILUT; \
  typedef Ifpack2::AdditiveSchwarz<MAT_P,PRECOND_ILUT> PRECOND_AS; \
  typedef Ifpack2::Diagonal<MAT_P> PRECOND_D; \
  typedef Tpetra::details::ApplyOp<Scalar,PRECOND_ILUT> PRECOND_ILUT_OP; \
  typedef Tpetra::details::ApplyOp<Scalar,PRECOND_AS> PRECOND_AS_OP;	\
  typedef Tpetra::details::ApplyOp<Scalar,PRECOND_D> PRECOND_D_OP;

#elif MIXED_PREC == 0

#define TYPEDEF_MIXED(SCALAR, LO, GO, NODE) \
  \
  typedef Scalar precScalar; \
  typedef typename std::map<GO, precScalar>::iterator ITER; \
  typedef MAT MAT_P; \
  typedef MMOP MMOP_P; \
  typedef OP OP_P; \
  typedef VEC VEC_P; \
  typedef MV MV_P; \
  typedef PRECOND PRECOND_P; \
  typedef SCALE SCALE_P; \
  typedef Ifpack2::ILUT<MAT_P> PRECOND_ILUT; \
  typedef Ifpack2::AdditiveSchwarz<MAT_P,PRECOND_ILUT> PRECOND_AS; \
  typedef Ifpack2::Diagonal<MAT_P> PRECOND_D; \
  typedef Tpetra::details::ApplyOp<Scalar,PRECOND_ILUT> PRECOND_ILUT_OP; \
  typedef Tpetra::details::ApplyOp<Scalar,PRECOND_AS> PRECOND_AS_OP;	\
  typedef Tpetra::details::ApplyOp<Scalar,PRECOND_D> PRECOND_D_OP;

#endif // MIXED_PREC

#endif //SUPPORTS_STRATIMIKOS

#endif /* TPETRA_HEADERS_H */
