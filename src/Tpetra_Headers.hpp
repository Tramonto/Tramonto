#ifndef TPETRA_HEADERS_H
#define TPETRA_HEADERS_H

#include "Kokkos_NodeTrace.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_OperatorApplyInverse.hpp"
#include "Tpetra_InvOperator.hpp"
#include "Tpetra_ParameterListConverter.hpp"
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
#include <MueLu.hpp>
#include <MueLu_TrilinosSmoother.hpp>
#include "MueLu_Hierarchy.hpp"
#include <MueLu_RAPFactory.hpp>
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Amesos2Smoother.hpp"
#include "MueLu_Utilities.hpp"
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>
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
  typedef Teuchos::Comm<int> COMM; \
  typedef Tpetra::Map<LO,GO,Node> MAP; \
  typedef Tpetra::CrsMatrix<SCALAR,LO,GO,Node> MAT; \
  typedef Tpetra::Import<LO,GO,Node> IMP; \
  typedef Belos::SolverManager<SCALAR, MV, OP> SolMGR; \
  typedef Belos::LinearProblem<SCALAR, MV, OP> LinPROB; \
  typedef Ifpack2::Preconditioner<SCALAR, LO, GO, Node> PRECOND; \
  typedef typename Teuchos::ScalarTraits<SCALAR>::halfPrecision halfScalar; \
  typedef Tpetra::MultiVector<halfScalar,LO,GO,Node> MV_H; \
  typedef Tpetra::Vector<halfScalar,LO,GO,Node> VEC_H;		   \
  typedef Tpetra::Operator<halfScalar,LO,GO,Node> OP_H;			\
  typedef Tpetra::OperatorApplyInverse<halfScalar,LO,GO,Node> APINV_H;	\
  typedef Tpetra::InvOperator<halfScalar,LO,GO,Node> INVOP_H; \
  typedef Tpetra::CrsMatrix<halfScalar,LO,GO,Node> MAT_H;			\
  typedef Belos::SolverManager<halfScalar, MV_H, OP_H> SolMGR_H;		\
  typedef Belos::LinearProblem<halfScalar, MV_H, OP_H> LinPROB_H;		\
  typedef Ifpack2::Preconditioner<halfScalar, LO, GO, Node> PRECOND_H;

#endif //SUPPORTS_STRATIMIKOS

#endif /* TPETRA_HEADERS_H */
