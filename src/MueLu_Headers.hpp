#if ENABLE_MUELU == 1
#include <MueLu_TrilinosSmoother.hpp>
#include <MueLu_Hierarchy.hpp>
#include <MueLu_Amesos2Smoother.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_MLParameterListInterpreter.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_UseShortNames.hpp>
#endif

#define MUELU_TYPEDEF(MSC, LO, GO, NO, LMO)	\
  \
  typedef MueLu::Hierarchy<MSC,LO,GO,NO,LMO> Hierarchy; \
  typedef MueLu::MLParameterListInterpreter<MSC,LO,GO,NO,LMO> MLFactory; \
  typedef Xpetra::TpetraCrsMatrix<MSC,LO,GO,NO,LMO> XpetraTpetraCrsMatrix; \
  typedef Xpetra::CrsMatrix<MSC,LO,GO,NO,LMO> XpetraCrsMatrix; \
  typedef Xpetra::Matrix<MSC,LO,GO,NO,LMO> XpetraMatrix; \
  typedef Xpetra::CrsMatrixWrap<MSC,LO,GO,NO,LMO> XpetraCrsWrap; \
  typedef MueLu::TpetraOperator<MSC,LO,GO,NO,LMO> MueLuOP;
