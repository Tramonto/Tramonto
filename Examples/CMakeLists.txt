#
# These are the primary Tramonto tests.
#

include(CMakeParseArguments)

# Set a cache variable to serve as the default for the maximum number
# of processes a test can use.
set(Tramonto_TEST_HARNESS_DEFAULT_PROCESSES 4 CACHE STRING 
  "Maximum number of processes to be used by individual test harness jobs.")

# We use MPI to run Tramonto tests.
find_package(MPI REQUIRED)

# Add a Tramonto executable test.
#
# This function adds a test that runs the Tramonto executable with a
# set of input files. The input files are specified with the
# INPUT_FILES name argument. A set of regular expressions to test for
# in the program output is specified in the README file.
function(Tramonto_add_test _test_name)

  # Set up the accepted argument set for the macro.
  set(_input_options "")
  set(_input_one_value_keywords README_FILE PROCESSES)
  set(_input_multi_value_keywords INPUT_FILES)

  # Parse the arguments.
  cmake_parse_arguments(
    TEST_HARNESS
    ${_input_options} 
    ${_input_one_value_keywords}
    ${_input_multi_value_keywords}
    ${ARGN})

  
  # Set defaults for required values.
  if(NOT TEST_HARNESS_README_FILE)
    set(TEST_HARNESS_README_FILE "README")
  endif()

  if(NOT TEST_HARNESS_PROCESSES)
    set(TEST_HARNESS_PROCESSES "${Tramonto_TEST_HARNESS_DEFAULT_PROCESSES}")
  endif()

  # Add targets to copy all of the test dependencies to the binary
  # directory.
  foreach(_file ${TEST_HARNESS_INPUT_FILES})
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_file}
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_CURRENT_SOURCE_DIR}/${_file} 
	${CMAKE_CURRENT_BINARY_DIR}/${_file}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    add_custom_target(copy_${_test_name}_${_file} ALL
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()

  # Add the test.
  get_property(_dft_exec_location TARGET dft_exec PROPERTY LOCATION)
  add_test(
    NAME ${_test_name}
    COMMAND 
      ${MPIEXEC} 
      ${MPIEXEC_NUMPROC_FLAG} ${TEST_HARNESS_PROCESSES}
      ${MPIEXEC_PREFLAGS} ${_dft_exec_location} ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  
  # Set properties on the test.
  set_tests_properties(
    ${_test_name} 
    PROPERTIES
      PROCESSORS ${TEST_HARNESS_PROCESSES}
      WILL_FAIL TRUE
      )
  
endfunction()

add_subdirectory(HS1_FMT1_1D)