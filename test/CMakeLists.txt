add_subdirectory(third_party)

# add the unit tests in the subdirectories
add_test_files( TEST_FILES "ut" lib/adaptation )
add_test_files( TEST_FILES "ut" lib/common )
add_test_files( TEST_FILES "ut" lib/geometry )
add_test_files( TEST_FILES "ut" lib/graphics )
add_test_files( TEST_FILES "ut" lib/library )
add_test_files( TEST_FILES "ut" lib/mesh )
add_test_files( TEST_FILES "ut" lib/mesh/delaunay )
add_test_files( TEST_FILES "ut" lib/mesh/frontal )
add_test_files( TEST_FILES "ut" lib/numerics )
add_test_files( TEST_FILES "ut" lib/shape )
add_test_files( TEST_FILES "ut" lib/voronoi )
add_test_files( TEST_FILES "ut" bin )

# sandbox files should have _toy.cpp
add_test_files( SBX_FILES "toy" sandbox/philip )
add_test_files( SBX_FILES "toy" sandbox/hannah )
add_test_files( SBX_FILES "toy" sandbox/paul )

set( unit_skip mesh/delaunay/vertex_ut.cpp )

# the main unit test file
set( MASTER_UNIT ut_driver.cpp )
set( COMMON_UNIT ut_common.cpp )

add_library( ut_common_lib STATIC ${COMMON_UNIT} )
target_link_libraries( ut_common_lib avro_shared )

# needed for the directory class definitions
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/third_party)
include_directories( ${CMAKE_SOURCE_DIR}/src/third_party/numpack )
include_directories( third_party/glad/include )
include_directories( third_party/glfw/include )

set( CHECK_TESTS )

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
  target_link_libraries( ${test_exe} ut_common_lib )
  target_compile_definitions( ${test_exe} PUBLIC -DSTANDALONE )

  # TODO re-add parallel testing capability...but how...
  # target that actually runs the executable
  #list( GET TEST_NUMPROCS ${i} numprocs )

  #if ( ${numprocs} GREATER 0 AND avro_RUN_WITH_MPI)
#    RUN_WITH_MPIEXEC( ${numprocs} )
#  else()
    set( MPI_COMMAND )
#  endif()
  add_custom_target( ${test_bin} command ${MPI_COMMAND} $<TARGET_FILE:${test_exe}> 1 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )

  add_custom_target( ${test_bin}_off command ${MPI_COMMAND} $<TARGET_FILE:${test_exe}> 0 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )
  set( check_test ${test_bin}_off )

  set( CHECK_TESTS ${CHECK_TESTS} ${check_test} )

  # target for coverage: cmake/coverage.cmake
  ADD_COVERAGE_UT( ${test_bin}_coverage ${test_bin} )

  # targets for memory & optionally stack checking & gdb
  ADD_DEBUG_TARGETS( ${test_bin} ${CMAKE_SOURCE_DIR}/test )

endforeach()

# create the full unit executable
avro_add_executable( unit_exe ${MASTER_UNIT} ${TEST_FILES} )
target_link_libraries( unit_exe avro_shared )
target_compile_definitions( unit_exe PUBLIC -DSTDOUT_REDIRECT="${CMAKE_BINARY_DIR}/unit_tests_output.txt" )

# target for running all unit tests independently
add_custom_target( check_tests )
add_dependencies( check_tests ${CHECK_TESTS} )

# target for running all unit tests from a single executable
add_custom_target( unit command $<TARGET_FILE:unit_exe> WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )
ADD_DEBUG_TARGETS( unit ${CMAKE_SOURCE_DIR}/test )
ADD_COVERAGE_UT( unit_coverage unit )

# create targets for the sandbox files
set( SBX_COMPILE )
foreach( test_src ${SBX_FILES} )

  # replace the '/' with '_'
  string( REPLACE "/" "_" test_src_underscored ${test_src} )

  # remove the extension
  get_filename_component( test_bin ${test_src_underscored} NAME_WE )

  # assign the name of the executable
  set( test_exe ${test_bin}_exe )

  # the actual executable: target for simply building
  avro_add_executable( ${test_exe} ${test_src} )
  target_link_libraries( ${test_exe} ut_common_lib )
  target_compile_definitions( ${test_exe} PUBLIC -DSTANDALONE )

  # add the target to all sandbox targets
  set( SBX_COMPILE ${SBX_COMPILE} ${test_exe} )

  # target for running the sandbox case with the plotter turned on
  add_custom_target( ${test_bin} command $<TARGET_FILE:${test_exe}> 1 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )

  # targets for memory & optionally stack checking & gdb
  ADD_DEBUG_TARGETS( ${test_bin} ${CMAKE_SOURCE_DIR}/test )

endforeach()

# target for compiling all sandbox tests
add_custom_target( sandbox )
add_dependencies( sandbox ${SBX_COMPILE} )

# add a target for the regression tests
add_custom_target( check_benchmarks ${CMAKE_SOURCE_DIR}/test/regression/benchmark3d.sh WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test )

# reset the directory for the executables
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
