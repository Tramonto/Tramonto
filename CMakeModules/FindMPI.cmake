# Override of FindMPI to make use of MPI settings from Trilinos. If
# Trilinos cannot be found, we default back to the system version of
# FindMPI.

find_package(Trilinos REQUIRED PATHS ${TRILINOS_PATH}/include ${TRILINOS_PATH})
if(Trilinos_FOUND)
  # Set the FindMPI variables with their Trilinos equivalents.
  set(MPI_C_INCLUDE_PATH "${Trilinos_MPI_INCLUDE_DIRS}" CACHE STRING "From Trilinos")
  set(MPI_C_LIBRARIES "${Trilinos_MPI_LIBRARIES}" CACHE STRING "From Trilinos")
  set(MPIEXEC "${Trilinos_MPI_EXEC}" CACHE STRING "From Trilinos")
  set(MPIEXEC_NUMPROC_FLAG "${Trilinos_MPI_EXEC_NUMPROCS_FLAG}" CACHE STRING "From Trilinos")
  set(MPIEXEC_PREFLAGS "" CACHE STRING "From Trilinos")
  set(MPIEXEC_POSTFLAGS "" CACHE STRING "From Trilinos")
endif()

include(${CMAKE_ROOT}/Modules/FindMPI.cmake)
