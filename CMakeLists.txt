cmake_minimum_required(VERSION 2.6)

MESSAGE("Starting configure")
FIND_PACKAGE(Trilinos PATHS ${TRILINOS_PATH})

#Verify that we can find a Trilinos installation
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
  MESSAGE("Could not find Trilinos 10.2 or greater!")
ENDIF()

IF(NOT CMAKE_CXX_COMPILER)
  MESSAGE("Setting CXX compiler.")
  SET(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
ENDIF()

IF(NOT CMAKE_C_COMPILER)
  MESSAGE("Setting C compiler.")
  SET(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
ENDIF()

IF(NOT CMAKE_Fortran_COMPILER)
  MESSAGE("Setting Fortran compiler.")
  SET(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
ENDIF()

IF(NOT CMAKE_CXX_COMPILER_FLAGS)
  MESSAGE("Setting CXX compiler flags.")
  SET(CMAKE_CXX_COMPILER_FLAGS ${Trilinos_CXX_COMPILER_FLAGS})
ENDIF()

IF(NOT CMAKE_C_COMPILER_FLAGS)
  MESSAGE("Setting C compiler flags.")
  SET(CMAKE_C_COMPILER_FLAGS ${Trilinos_C_COMPILER_FLAGS})
ENDIF()

IF(NOT CMAKE_Fortran_COMPILER_FLAGS)
  MESSAGE("Setting Fortran compiler flags.")
  SET(CMAKE_Fortran_COMPILER_FLAGS ${Trilinos_Fortran_COMPILER_FLAGS})
ENDIF()

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
#enable_language(Fortran)
#Currently no example exes are built for Tramonto
#set(BUILD_EXAMPLES true CACHE BOOL "Build Tramonto Examples")
set(TRILINOS_PATH "" CACHE FILEPATH "Path to Trilinos Directory")
#include(FindMPI)
#include(FindBLAS)
#include(FindLAPACK)

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

