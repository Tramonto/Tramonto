cmake_minimum_required(VERSION 2.6)
project(Tramonto)
enable_language(Fortran)
#Currently no example exes are built for Tramonto
#set(BUILD_EXAMPLES true CACHE BOOL "Build Tramonto Examples")
set(TRILINOS_PATH "" CACHE FILEPATH "Path to Trilinos Directory")
include(FindMPI)
include(FindBLAS)
include(FindLAPACK)

#if(BUILD_EXAMPLES)
#	ADD_SUBDIRECTORY(Examples)
#endif(BUILD_EXAMPLES)
if(NOT TRILINOS_PATH)
	message(FATAL_ERROR "You must specifiy the path to Trilinos")
endif(NOT TRILINOS_PATH)

if(NOT MPI_FOUND)
	message(FATAL_ERROR "Could not find mpi")
endif(NOT MPI_FOUND)
if(NOT BLAS_FOUND)
	message(FATAL_ERROR "Could not find blas")
endif(NOT BLAS_FOUND)
if(NOT LAPACK_FOUND)
	message(FATAL_ERROR "Could not find lapack")
endif(NOT LAPACK_FOUND)

	
ADD_SUBDIRECTORY(src)

