cmake_minimum_required(VERSION 2.6)

MESSAGE("Starting configure")

# Use TRILINOS_PATH if the user set it, to help find Trilinos 10.7 or greater.
# The final location will actually be held in Trilinos_DIR which must
# point at "<prefix>/lib/cmake/Trilinos", but this helps the search.
SET(CMAKE_PREFIX_PATH ${TRILINOS_PATH} ${CMAKE_PREFIX_PATH})

#Look for Trilinos - the first PATH will find Trilinos 10.4 and 10.6,
#the second PATH will find Trilinos 10.2
FIND_PACKAGE(Trilinos REQUIRED PATHS ${TRILINOS_PATH}/include ${TRILINOS_PATH})

#Output some info about the Trilinos installation found
IF(Trilinos_FOUND)
   MESSAGE("\nFound Trilinos!  Here are the details: ")
   MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
   MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
   MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
   MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
   MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
   MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
   MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
   # Will Tramonto work with a Shared version of Trilinos?
   MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
   #These last outputs are temporary
   MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
   MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
   MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
   MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS =
     ${Trilinos_Fortran_COMPILER_FLAGS}")
   MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
   MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
   MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
   # Add MPI_EXEC, MPI_EXEC_MAX_NUMPROCS, MPI_EXEC_NUMPROC_FLAGS
   # These are not available in Trilinos 10.2.0, but should be in 10.4.0
   MESSAGE("End of Trilinos details\n")
ELSE()
  # We should never reach this because FIND_PACKAGE(Trilinos now uses REQUIRED
  MESSAGE(FATAL_ERROR "Could not find Trilinos 10.2 or greater!")
ENDIF()

IF(NOT CMAKE_CXX_COMPILER)
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
ENDIF()

IF(NOT CMAKE_C_COMPILER)
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
ENDIF()

#IF(NOT CMAKE_Fortran_COMPILER)
#  MESSAGE("Setting Fortran compiler.")
#  SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
#ENDIF()

IF(NOT CMAKE_LINKER)
  MESSAGE("Setting linker.")
  SET(CMAKE_LINKER ${Trilinos_LINKER})
ENDIF()

IF(NOT CMAKE_EXE_LINKER_FLAGS)
  MESSAGE("Setting linker flags.")
  SET(CMAKE_EXE_LINKER_FLAGS ${Trilinos_EXTRA_LD_FLAGS})
ENDIF()

IF(NOT CMAKE_AR)
  MESSAGE("Setting ar.")
  SET(CMAKE_AR ${CMAKE_AR})
ENDIF()

project(Tramonto)
SET(CPACK_PACKAGE_VERSION "4.0")
#enable_language(Fortran)
#Currently no example exes are built for Tramonto
#set(BUILD_EXAMPLES true CACHE BOOL "Build Tramonto Examples")
set(TRILINOS_PATH "" CACHE FILEPATH "Path to Trilinos Directory")
#include(FindMPI)
#include(FindBLAS)
#include(FindLAPACK)

IF(NOT CMAKE_CXX_FLAGS)
  MESSAGE("Setting CXX compiler flags.")
  SET(CMAKE_CXX_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})
ENDIF()

IF(NOT CMAKE_C_FLAGS)
  MESSAGE("Setting C compiler flags.")
  SET(CMAKE_C_FLAGS ${Trilinos_C_COMPILER_FLAGS})
ENDIF()

#IF(NOT CMAKE_Fortran_FLAGS)
#  MESSAGE("Setting Fortran compiler flags.")
#  SET(CMAKE_Fortran_FLAGS ${Trilinos_Fortran_COMPILER_FLAGS})
#ENDIF()

#if(BUILD_EXAMPLES)
#	ADD_SUBDIRECTORY(Examples)
#endif(BUILD_EXAMPLES)
#if(NOT TRILINOS_PATH)
#	message(FATAL_ERROR "You must specifiy the path to Trilinos")
#endif(NOT TRILINOS_PATH)

#if(NOT MPI_FOUND)
#	message(FATAL_ERROR "Could not find mpi")
#endif(NOT MPI_FOUND)
#if(NOT BLAS_FOUND)
#	message(FATAL_ERROR "Could not find blas")
#endif(NOT BLAS_FOUND)
#if(NOT LAPACK_FOUND)
#	message(FATAL_ERROR "Could not find lapack")
#endif(NOT LAPACK_FOUND)

	
ADD_SUBDIRECTORY(src)

install(FILES
	Examples/BC1_FMT1_WCA_ELEC_2D/dft_input.dat
	Examples/BC1_FMT1_WCA_ELEC_2D/dft_surfaces.dat
	Examples/BC1_FMT1_WCA_ELEC_2D/README
	DESTINATION ./Examples/BC1_FMT1_WCA_ELEC_2D)
Install(FILES
	Examples/BULK1_FMT1_WCA_1D/dft_input.dat
	Examples/BULK1_FMT1_WCA_1D/README
	DESTINATION ./Examples/BULK1_FMT1_WCA_1D)
Install(FILES
	Examples/CONT1_FMT1_WCA_1D/dft_input.dat
	Examples/CONT1_FMT1_WCA_1D/dft_surfaces.dat
	Examples/CONT1_FMT1_WCA_1D/README
	DESTINATION ./Examples/CONT1_FMT1_WCA_1D)
Install(FILES
	Examples/HS1_FMT1_1D/dft_input.dat
	Examples/HS1_FMT1_1D/dft_surfaces.dat
	Examples/HS1_FMT1_1D/README
	DESTINATION ./Examples/HS1_FMT1_1D)
Install(FILES
	Examples/LJ1_FMT1_WCA_1D/dft_input.dat
	Examples/LJ1_FMT1_WCA_1D/dft_surfaces.dat
	Examples/LJ1_FMT1_WCA_1D/README
	DESTINATION ./Examples/LJ1_FMT1_WCA_1D)
Install(FILES
	Examples/MIX1_FMT1_1D/dft_input.dat
	Examples/MIX1_FMT1_1D/dft_surfaces.dat
	Examples/MIX1_FMT1_1D/README
	DESTINATION ./Examples/MIX1_FMT1_1D)
Install(FILES
	Examples/NLSOLVER1_WJDC_FMT3_1D/dft_input.dat
	Examples/NLSOLVER1_WJDC_FMT3_1D/dft_surfaces.dat
	Examples/NLSOLVER1_WJDC_FMT3_1D/lin_18
	Examples/NLSOLVER1_WJDC_FMT3_1D/README
	DESTINATION ./Examples/NLSOLVER1_WJDC_FMT3_1D)
Install(FILES
	Examples/POLY1_CMS_1D/dft_input.dat
	Examples/POLY1_CMS_1D/dft_surfaces.dat
	Examples/POLY1_CMS_1D/README
	Examples/POLY1_CMS_1D/lin_18
	Examples/POLY1_CMS_1D/cr8.2.8_0.711eff1
	DESTINATION ./Examples/POLY1_CMS_1D)
Install(FILES
	Examples/POLYBRANCH1_WJDC_FMT3_1D/dft_input.dat
	Examples/POLYBRANCH1_WJDC_FMT3_1D/dft_surfaces.dat
	Examples/POLYBRANCH1_WJDC_FMT3_1D/README
	Examples/POLYBRANCH1_WJDC_FMT3_1D/star_3arm_5mer
	Examples/POLYBRANCH1_WJDC_FMT3_1D/dens.eta1
	Examples/POLYBRANCH1_WJDC_FMT3_1D/dens.eta2
	Examples/POLYBRANCH1_WJDC_FMT3_1D/dens.eta3
	DESTINATION ./Examples/POLYBRANCH1_WJDC_FMT3_1D)
Install(FILES
	Examples/SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D/dft_input.dat
	Examples/SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D/dft_dens.dat
	Examples/SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D/dft_input.dat
	Examples/SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D/README
	Examples/SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D/lin_18_1
	Examples/SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D/dens_pure_solvent.dat
	DESTINATION ./Examples/SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D)
Install(FILES
	Examples/YUKAWA_FMT1_3D/dft_input.dat
	Examples/YUKAWA_FMT1_3D/dft_surfaces.dat
	Examples/YUKAWA_FMT1_3D/README
	DESTINATION ./Examples/YUKAWA_FMT1_3D)
include(CPack)
