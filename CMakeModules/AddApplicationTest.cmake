# Add an application level test.
#
# An application level test consists of a preprocessing step, running
# the application itself, and a postprocessing step. In addition, file
# dependencies can be specified and an overall job weight can be
# specified.

include(CMakeParseArguments)

get_filename_component(_add_application_test_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

function(add_application_test _test_name)
  set(_option_args)
  set(_one_value_keyword_args 
    WEIGHT
    TIMEOUT
    WORKING_DIRECTORY)
  set(_multi_value_keyword_args
    FILE_DEPENDENCIES
    PREPROCESS
    APPLICATION 
    POSTPROCESS
    LABELS)

  cmake_parse_arguments(
    _application_test
    "${_option_args}"
    "${_one_value_keyword_args}"
    "${_multi_value_keyword_args}"
    ${ARGN})

  # Setup default values for the arguments.
  if(NOT _application_test_WEIGHT)
    set(_application_test_WEIGHT 1)
  endif()
  if(NOT _application_test_TIMEOUT)
    set(_application_test_TIMEOUT 0)
  endif()
  if(NOT _application_test_WORKING_DIRECTORY)
    set(_application_test_WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  # Add targets to copy all of the test dependencies to the working
  # directory at build time.
  foreach(_file ${_application_test_FILE_DEPENDENCIES})
    get_filename_component(_filename ${_file} NAME)
    
    add_custom_command(
      OUTPUT ${_application_test_WORKING_DIRECTORY}/${_file}
      COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_CURRENT_SOURCE_DIR}/${_file} 
	${_application_test_WORKING_DIRECTORY}/${_file}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    
    add_custom_target(copy_${_test_name}_${_filename} ALL
      DEPENDS ${_application_test_WORKING_DIRECTORY}/${_file})
  endforeach()

  # Parse each of the PREPROCESS, APPLICATION, POSTPROCESS argument lists.
  set(_phases PREPROCESS APPLICATION POSTPROCESS)
  foreach(_phase ${_phases})
    cmake_parse_arguments(
      _application_test_${_phase}
      ""
      "EXIT_CODE"
      "COMMAND"
      "${_application_test_${_phase}}")
    if(NOT _application_test_${_phase}_EXIT_CODE)
      set(_application_test_${_phase}_EXIT_CODE 0)
    endif()
  endforeach()

  # Configure the driver script.
  set(_test_driver_file ${_application_test_WORKING_DIRECTORY}/TEST-${_test_name}.cmake)
  configure_file(
    ${_add_application_test_dir}/ApplicationTest.cmake.in
    ${_test_driver_file}
    @ONLY)
  
  # Add the test that runs the driver.
  add_test(NAME ${_test_name}
    WORKING_DIRECTORY ${_application_test_WORKING_DIRECTORY}
    COMMAND ${CMAKE_COMMAND} -P ${_test_driver_file})
  
  # Set up the test properties.
  set_tests_properties( 
    ${_test_name}
    PROPERTIES
      PROCESSORS ${_application_test_WEIGHT}
      TIMEOUT ${_application_test_TIMEOUT}
      LABELS ${_application_test_WEIGHT})
    

endfunction()