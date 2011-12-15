#
# These are the primary Tramonto tests.
#

include(CMakeParseArguments)

# Set a cache variable to serve as the default for the maximum number
# of processes a test can use.
set(Tramonto_TEST_HARNESS_DEFAULT_PROCESSES 4 CACHE STRING 
  "Maximum number of processes to be used by individual test harness jobs.")

set(Tramonto_TEST_HARNESS_FAIL_TOLERANCE 0.001 CACHE STRING
  "Maximum difference between floating point values before test is considered a failure.")

# We use MPI to run Tramonto tests.
find_package(MPI REQUIRED)

# Build the program that checks the results.
add_subdirectory(check_results)

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
    "${_input_options}" 
    "${_input_one_value_keywords}"
    "${_input_multi_value_keywords}"
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
  set(_files_to_copy)
  list(APPEND  _files_to_copy ${TEST_HARNESS_INPUT_FILES} ${TEST_HARNESS_README_FILE})
  foreach(_file ${_files_to_copy})
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_file}
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_CURRENT_SOURCE_DIR}/${_file} 
	${CMAKE_CURRENT_BINARY_DIR}/${_file}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    add_custom_target(copy_${_test_name}_${_file} ALL
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()

  # Add the test that runs the example.
  get_property(_dft_exec_location TARGET dft_exec PROPERTY LOCATION)
  add_test(
    NAME ${_test_name}.run
    COMMAND 
      ${MPIEXEC} 
      ${MPIEXEC_NUMPROC_FLAG} ${TEST_HARNESS_PROCESSES}
      ${MPIEXEC_PREFLAGS} ${_dft_exec_location} ${MPIEXEC_POSTFLAGS}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  # Add the test that checks the results.
  add_test(
    NAME ${_test_name}.check
    COMMAND check_results README dft_output.dat ${Tramonto_TEST_HARNESS_FAIL_TOLERANCE}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  
  # Set properties on the test.
  set_tests_properties(
    ${_test_name}.run 
    PROPERTIES
      PROCESSORS ${TEST_HARNESS_PROCESSES}
      WILL_FAIL TRUE
      )
    
  # Set up the dependency between the run job and the check job.
  set_property(TEST ${_test_name}.check 
    APPEND PROPERTY DEPENDS ${_test_name}.run)
  
endfunction()

add_subdirectory(HS1_FMT1_1D)
