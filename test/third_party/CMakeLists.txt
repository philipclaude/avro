include_directories( ${CMAKE_SOURCE_DIR}/test )

# add the unit tests in the subdirectories
add_test_files( TEST_FILES "ut" numpack/DenseLinAlg )

# the master unit test file
set( MASTER_UNIT ${CMAKE_SOURCE_DIR}/test/ut_driver.cpp )
set( COMMON_UNIT ${CMAKE_SOURCE_DIR}/test/ut_common.cpp )

add_library( ut_third_party_common_lib STATIC ${COMMON_UNIT} )
target_link_libraries( ut_third_party_common_lib avro_shared )

# needed for the directory class definitions
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories( ${CMAKE_SOURCE_DIR}/src/third_party/numpack )

# set the directory for the executables
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/test/bin )

# create a target for each test source
foreach( test_src ${TEST_FILES} )

  list (FIND unit_skip ${test_src} _index)

  # replace the '/' with '_'
  string( REPLACE "/" "_" test_src_underscored ${test_src} )

  # remove the extension
  get_filename_component( test_bin ${test_src_underscored} NAME_WE )

  # assign the name of the executable
  set( test_exe ${test_bin}_exe )

  # the actual executable: target for simply building
  avro_add_executable( ${test_exe} ${test_src} )
  target_link_libraries( ${test_exe} ut_third_party_common_lib )
  target_compile_definitions( ${test_exe} PUBLIC -DSTANDALONE )

  add_custom_target( ${test_bin} command ${MPI_COMMAND} $<TARGET_FILE:${test_exe}> 1 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )

  add_custom_target( ${test_bin}_off command ${MPI_COMMAND} $<TARGET_FILE:${test_exe}> 0 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )
  set( check_test ${test_bin}_off )

  # target for coverage: cmake/coverage.cmake
  ADD_COVERAGE_UT( ${test_bin}_coverage ${test_bin} )

  # targets for memory & optionally stack checking & gdb
  ADD_DEBUG_TARGETS( ${test_bin} ${CMAKE_SOURCE_DIR}/test )

endforeach()

# create the full unit executable
avro_add_executable( third_party_unit_exe ${MASTER_UNIT} ${TEST_FILES} )
target_link_libraries( third_party_unit_exe avro_shared )
target_compile_definitions( third_party_unit_exe PUBLIC -DSTDOUT_REDIRECT="${CMAKE_BINARY_DIR}/unit_tests.out" )

# target for running all unit tests from a single executable
add_custom_target( unit_third_party command $<TARGET_FILE:unit_exe> WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )

# reset the directory for the executables
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
