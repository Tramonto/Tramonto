#
# These are the primary Tramonto tests.
#

include(AddApplicationTest)
include(CMakeParseArguments)
include(ConfigureTestOutputSize)

# Set a cache variable to serve as the default for the maximum number
# of processes a test can use.
set(Tramonto_TEST_HARNESS_DEFAULT_PROCESSES 4 CACHE STRING 
  "Maximum number of processes to be used by individual test harness jobs.")

set(Tramonto_TEST_HARNESS_FAIL_TOLERANCE 0.001 CACHE STRING
  "Maximum difference between floating point values before test is considered a failure.")

set(Tramonto_TEST_HARNESS_ITERATION_BOUND 4 CACHE STRING
  "Maximum number of iterations a test is allowed to differ by before it is considered a failure.")

set(Tramonto_TEST_HARNESS_Sanity_TIMEOUT 20 CACHE STRING
  "Maximum runtime for a test that gets run as part of the Sanity test suite.")

set(Tramonto_TEST_HARNESS_Checkin_TIMEOUT 30 CACHE STRING
  "Maximum runtime for a test that gets run as part of the Checkin test suite.")

set(Tramonto_TEST_HARNESS_Continuous_TIMEOUT 300 CACHE STRING
  "Maximum runtime for a test that gets run as part of the Continuous test suite.")

set(Tramonto_TEST_HARNESS_Nightly_TIMEOUT 3600 CACHE STRING
  "Maximum runtime for a test that gets run as part of the Nightly test suite.")

set(Tramonto_TEST_HARNESS_Occasional_TIMEOUT 36000 CACHE STRING
  "Maximum runtime for a test that gets run as part of the Occasional test suite.")

set(Tramonto_TEST_HARNESS_CHECK_JOB_TIMEOUT 2 CACHE STRING
  "Maximum runtime for the jobs responsible for checking the results of the test harness.")

# We use MPI to run Tramonto tests.
find_package(MPI REQUIRED)

# Build the program that checks the results.
add_subdirectory(check_results)

# Set the test output size to unlimited.
configure_test_output_size(0 0)

# Add a Tramonto executable test.
#
# This function adds a test that runs the Tramonto executable with a
# set of input files. The input files are specified with the
# INPUT_FILES name argument. A set of regular expressions to test for
# in the program output is specified in the README file.
function(Tramonto_add_test _test_name)

  # Set up the accepted argument set for the macro.
  set(_input_options)
  set(_input_one_value_keywords PROCESSES TYPE)
  set(_input_multi_value_keywords INPUT_FILES)

  # Parse the arguments.
  cmake_parse_arguments(
    TEST_HARNESS
    "${_input_options}" 
    "${_input_one_value_keywords}"
    "${_input_multi_value_keywords}"
    ${ARGN})

  
  # Set defaults for required values.
  if(NOT TEST_HARNESS_TYPE)
    set(TEST_HARNESS_TYPE Checkin)
  endif()

  if(NOT TEST_HARNESS_PROCESSES)
    set(TEST_HARNESS_PROCESSES "${Tramonto_TEST_HARNESS_DEFAULT_PROCESSES}")
  endif()

  if(NOT TEST_HARNESS_INPUT_FILES)
    file(GLOB TEST_HARNESS_INPUT_FILES *)
  endif()

  # Add the test that runs the example.
  get_property(_dft_exec_location TARGET dft_exec PROPERTY LOCATION)
  set(_command
    ${MPIEXEC} 
    ${MPIEXEC_NUMPROC_FLAG} ${TEST_HARNESS_PROCESSES}
    ${MPIEXEC_PREFLAGS} ${_dft_exec_location} ${MPIEXEC_POSTFLAGS})
  get_property(_check_results_location TARGET check_results PROPERTY LOCATION)
  set(_postprocess
    ${_check_results_location} README dft_output.dat ${Tramonto_TEST_HARNESS_FAIL_TOLERANCE} ${Tramonto_TEST_HARNESS_ITERATION_BOUND})
  
  # Set up test labels based on the specified testing track.
  set(_labels)
  set(_timeout ${Tramonto_TEST_HARNESS_Nightly_TIMEOUT})
  if(TEST_HARNESS_TYPE STREQUAL Sanity)
    set(_labels Sanity Checkin Continuous Nightly)
    set(_timeout ${Tramonto_TEST_HARNESS_Sanity_TIMEOUT})
  elseif(TEST_HARNESS_TYPE STREQUAL Checkin)
    set(_labels Checkin Continuous Nightly)
    set(_timeout ${Tramonto_TEST_HARNESS_Checkin_TIMEOUT})
  elseif(TEST_HARNESS_TYPE STREQUAL Continuous)
    set(_labels Continuous Nightly)
    set(_timeout ${Tramonto_TEST_HARNESS_Continuous_TIMEOUT})
  elseif(TEST_HARNESS_TYPE STREQUAL Nightly)
    set(_labels Nightly)
    set(_timeout ${Tramonto_TEST_HARNESS_Nightly_TIMEOUT})
  elseif(TEST_HARNESS_TYPE STREQUAL Occasional)
    set(_labels Occasional)
    set(_timeout ${Tramonto_TEST_HARNESS_Occasional_TIMEOUT})
  else()
    set(_labels ${TEST_HARNESS_TYPE})
  endif()
  
  # Add the test.
  add_application_test(${_test_name}
    WEIGHT ${TEST_HARNESS_PROCESSES}
    TIMEOUT ${_timeout}
    FILE_DEPENDENCIES ${TEST_HARNESS_INPUT_FILES}
    PREPROCESS COMMAND "\"${CMAKE_COMMAND}\"" -E remove -f dft_output.dat
    APPLICATION COMMAND ${_command} EXIT_CODE 1
    POSTPROCESS COMMAND ${_postprocess}
    LABELS ${_labels})

endfunction()

add_custom_target(sanity_test
  ${CMAKE_CTEST_COMMAND} -L Sanity)

# Add a target to run only the short tests before a checkin.
add_custom_target(checkin_test
  ${CMAKE_CTEST_COMMAND} -L Checkin
  )

add_custom_target(continuous_test
  ${CMAKE_CTEST_COMMAND} -L Continuous
  )

add_custom_target(nightly_test
  ${CMAKE_CTEST_COMMAND} -L Nightly
  )

# Begin example listing
add_subdirectory(BC10_FMT1_WCA_ELEC_3D)
add_subdirectory(BC11_FMT1_WCA_ELEC_3D)
add_subdirectory(BC12_FMT1_WCA_ELEC_3D)
add_subdirectory(BC1_FMT1_WCA_ELEC_2D)
add_subdirectory(BC2_FMT1_WCA_ELEC_2D)
add_subdirectory(BC3_FMT1_WCA_ELEC_2D)
add_subdirectory(BC4_FMT1_WCA_ELEC_2D)
add_subdirectory(BC5_FMT1_WCA_ELEC_2D)
add_subdirectory(BC6_FMT1_WCA_ELEC_2D)
add_subdirectory(BC7_FMT1_WCA_ELEC_3D)
add_subdirectory(BC8_FMT1_WCA_ELEC_3D)
add_subdirectory(BC9_FMT1_WCA_ELEC_3D)
add_subdirectory(BULK1_FMT1_WCA_1D)
add_subdirectory(BULK2_FMT1_WCA_1D)
add_subdirectory(BULK3_FMT1_WCA_1D)
add_subdirectory(BULK4_MIX_FMT3_WCA_WJDC_1D)
add_subdirectory(BULK5_MIX_FMT3_WCA_WJDC_1D)
add_subdirectory(BULK6_MIX_FMT3_WCA_WJDC_1D)
add_subdirectory(BULK7_FMT3_WCA_WJDC_1D)
add_subdirectory(BULK8_MIX_FMT3_WCA_WJDC_1D)
add_subdirectory(CONT1_FMT1_WCA_1D)
add_subdirectory(CONT2_FMT1_WCA_1D)
add_subdirectory(CONT3_FMT1_WCA_1D)
add_subdirectory(CONT4_FMT1_WCA_1D)
add_subdirectory(DIFFUSION1_FMT1_WCA_1D)
add_subdirectory(DIFFUSION2_FMT1_WCA_WJDC_1D)
add_subdirectory(DIFFUSION3_FMT1_WCA_PICARD_1D)
add_subdirectory(DIFFUSION4_FMT1_WCA_3D)
add_subdirectory(DIFFUSION5_FMT1_WCA_WJDC_3D)
add_subdirectory(DIFFUSION6_FMT1_WCA_PICARD_3D)
add_subdirectory(EXP_FMT3_1D)
add_subdirectory(GEOM10_SPHERE_3D)
add_subdirectory(GEOM11_CYL_2D)
add_subdirectory(GEOM12_CYL_2D)
add_subdirectory(GEOM13_CYL_ROUGH_2D)
add_subdirectory(GEOM14_CYL_WEDGECUTOUT_2D)
add_subdirectory(GEOM15_CYL_WEDGE2TYPES_2D)
add_subdirectory(GEOM16_CYL_WEDGE2TYPES_REFLECT_2D)
add_subdirectory(GEOM17_1ELEM_1D)
add_subdirectory(GEOM18_1ELEM_2D)
add_subdirectory(GEOM19_1ELEM_3D)
add_subdirectory(GEOM1_PLANE_ROUGH_2D)
add_subdirectory(GEOM20_CYL_3D)
add_subdirectory(GEOM21_CYL_ROUGH_3D)
add_subdirectory(GEOM22_CYL_WEDGE2TYPES_3D)
add_subdirectory(GEOM23_CYL_PERIODIC_3D)
add_subdirectory(GEOM24_ATOMS_3D)
add_subdirectory(GEOM25_CYL_PORE_2D)
add_subdirectory(GEOM26_CYL_PORE_ROUGH_2D)
add_subdirectory(GEOM27_CYL_PORE_WEDGE2TYPES_2D)
add_subdirectory(GEOM28_CYLPORE_WEDGE2TYPES_3D)
add_subdirectory(GEOM29_SLITPORE_2D)
add_subdirectory(GEOM2_PLANE_ROUGH_2D)
add_subdirectory(GEOM30_TAPERED_SLITPORE_2D)
add_subdirectory(GEOM31_TAPERED_CYLPORE_3D)
add_subdirectory(GEOM32_PLANE_PERIODIC_2D)
add_subdirectory(GEOM33_PLANE_PERIODIC2_2D)
add_subdirectory(GEOM34_PLANE_PERIODIC3_2D)
add_subdirectory(GEOM35_PLANE_PERIODIC_3D)
add_subdirectory(GEOM36_PLANE_NOTCH_2D)
add_subdirectory(GEOM37_WEDGE_LINEAR_2D)
add_subdirectory(GEOM38_WEDGE_3D)
add_subdirectory(GEOM39_AMFTIP_3D)
add_subdirectory(GEOM3_PLANE_ROUGH_3D)
add_subdirectory(GEOM40_AMFTIP_HEMISPHERE_3D)
add_subdirectory(GEOM4_PLANE_ROUGH_3D)
add_subdirectory(GEOM5_PLANE_ROUGH_3D)
add_subdirectory(GEOM6_PLANE_ROUGH_VEXTINT_2D)
add_subdirectory(GEOM7_PLANE_ROUGH_VEXTINT_3D)
add_subdirectory(GEOM8_BLOCK_2D)
add_subdirectory(GEOM9_BLOCK_3D)
add_subdirectory(HS1_FMT1_1D)
add_subdirectory(HS2_FMT1_1D)
add_subdirectory(HS3_FMT1_1D)
add_subdirectory(HS4_FMT1_1D)
add_subdirectory(HS5_FMT1_1D)
add_subdirectory(HS6_FMT1_1D)
add_subdirectory(HS7_FMT1_1D)
add_subdirectory(HS8_FMT4_1D)
add_subdirectory(LJ1_FMT1_WCA_1D)
add_subdirectory(LJ2_FMT1_WCA_1D)
add_subdirectory(LJ3_FMT1_WCA_1D)
add_subdirectory(LJ4_FMT1_WCA_1D)
add_subdirectory(LJ5_FMT1_WCA2_1D)
add_subdirectory(MIX10_FMT1_ATTELEC_1D)
add_subdirectory(MIX11_FMT1_WCA_1D)
add_subdirectory(MIX12_FMT1_ELEC_1D)
add_subdirectory(MIX13_ELEC_PICARD_1D)
add_subdirectory(MIX1_FMT1_1D)
add_subdirectory(MIX2_FMT1_WCA_1D)
add_subdirectory(MIX3_ELEC_1D)
add_subdirectory(MIX4_FMT1_ELEC_1D)
add_subdirectory(MIX5_FMT1_ELEC_1D)
add_subdirectory(MIX6_FMT1_ELEC_1D)
add_subdirectory(MIX7_FMT1_WCA_ELEC_1D)
add_subdirectory(MIX8_FMT1_WCA_ELEC_1D)
add_subdirectory(MIX9_FMT1_WCA_ELEC_1D)
add_subdirectory(NLSOLVER1_WJDC_FMT3_1D)
add_subdirectory(NLSOLVER2_WJDC_FMT3_1D)
add_subdirectory(NLSOLVER3_WJDC_FMT3_1D)
add_subdirectory(NLSOLVER4_WJDC_FMT3_1D)
add_subdirectory(NLSOLVER5_WJDC_FMT3_1D)
add_subdirectory(POLY10_WJDC_FMT1_1D)
add_subdirectory(POLY11_WJDC_FMT1_1D)
add_subdirectory(POLY12_WJDC_FMT1_1D)
add_subdirectory(POLY13_WJDC_FMT1_WCA_1D)
add_subdirectory(POLY14_WJDC_FMT1_WCA_1D)
add_subdirectory(POLY15_WJDC_FMT3_1D)
add_subdirectory(POLY16_WJDC_FMT3_1D)
add_subdirectory(POLY17_WJDC_FMT3_1D)
add_subdirectory(POLY18_WJDC_FMT3_1D)
add_subdirectory(POLY19_WJDC_FMT3_1D)
add_subdirectory(POLY1_CMS_1D)
add_subdirectory(POLY20_CMS_1D)
add_subdirectory(POLY21_WJDC_FMT3_1D)
add_subdirectory(POLY22_WJDC_FMT3_1D)
add_subdirectory(POLY25_WJDC_FMT3_ELEC)
add_subdirectory(POLY26_WJDC_FMT3_ELEC)
add_subdirectory(POLY27_WJDC_FMT3_ELEC)
add_subdirectory(POLY28_WJDC_FMT3_ELEC)
add_subdirectory(POLY2_CMS_1D)
add_subdirectory(POLY3_CMS_1D)
add_subdirectory(POLY4_CMS_1D)
add_subdirectory(POLY5_WTC_FMT1_1D)
add_subdirectory(POLY6_WTC_FMT1_1D)
add_subdirectory(POLY7_WTC_FMT1_1D)
add_subdirectory(POLY8_WTC_FMT1_ELEC_1D)
add_subdirectory(POLY9_WTC_FMT1_1D)
add_subdirectory(POLYBRANCH1_WJDC_FMT3_1D)
add_subdirectory(POLYGRAFT1_WJDC_JAINPAPER_1D)
add_subdirectory(POLYGRAFT2_WJDC_1D)
add_subdirectory(POLYGRAFT3_WJDC_SEMIGRAND_1D)
add_subdirectory(POLYGRAFT4_WJDC_MIDSEG_1D)
add_subdirectory(POLYGRAFT5_WJDC_ATT_MIXTURE_1D)
add_subdirectory(SELF_ASSEMBLE1_FMT3_WCA_WJDC_1D)
add_subdirectory(SELF_ASSEMBLE2_FMT3_WCA_WJDC_1D)
add_subdirectory(SELF_ASSEMBLE3_FMT3_WCA_WJDC_1D)
add_subdirectory(SELF_ASSEMBLE4_FMT3_WCA_WJDC_1D)
add_subdirectory(SELF_ASSEMBLE5_BILAYERS_4MERTAIL)
add_subdirectory(SELF_ASSEMBLE6_BILAYERS_8MERTAIL)
add_subdirectory(SELF_ASSEMBLE7_BILAYERS_16MERTAIL)
add_subdirectory(SW_FMT3_1D)
add_subdirectory(YUKAWA_FMT1_3D)
