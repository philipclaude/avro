cmake_minimum_required(VERSION 2.8.12)

# this allows setting the version number
if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif()

# PkgConfig is used by many Finds/*.cmake
find_package(PkgConfig)

# setup libraries and build type
include(cmake/setup.cmake)
include(cmake/compiler.cmake)
project(avro VERSION 2.1)

if (POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

# default options
option(avro_STATIC_LINK_CXX "statically link to c++ libraries (probably because you're using voyager)" false)
option(AVRO_WITH_MPI "compile for distributed computing support" false)
option(RUN_WITH_MPI "run with mpi. might be difficult to debug with mpirun" false)
option(avro_WITH_FPG "compile FPG needed to make new predicates" false)
option(avro_NO_ESP "ESP is not available, so we will use internal wv and egads-lite" false )
set(avro_TEST_DATA "data/" CACHE STRING "location for avro test files relative to test/" )
option(avro_WITH_GL "option to compile OpenGL-based visualizer (along with WebGL-based one)" true)
option(avro_LITE "option to compile a lite version of avro" false )
option(AVRO_HEADLESS_GRAPHICS "option to use headless graphics (true for unit testing)" false)

set(avro_CPU_THREAD_MANAGER "cpp11" CACHE STRING "type of cpu threading environment to use [openmp, cpp11, pthread, emp]")
set(avro_GPU_THREAD_MANAGER "none" CACHE STRING "type of gpu threading environment to use [opencl,cuda]")

# schedule dynamic doesn't seem to be supported by intel
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  set(avro_CPU_THREAD_MANAGER "cpp11")
endif()

# avro requires the c++11 standard
set(CMAKE_CXX_STANDARD 11)

# overwrite anything for a lite version
if (${avro_LITE})
  message( STATUS "enabled lite version of avro" )
  set( AVRO_NO_ESP true )
  set( AVRO_WITH_MPI false )
  set( avro_WITH_GL false )
endif()

# required dependencies
include(cmake/FindLAPACK.cmake)
if (NOT AVRO_NO_ESP)
  include(cmake/FindESP.cmake)
endif()
include(cmake/FindNLOPT.cmake)
include(cmake/FindGraphics.cmake)

set(EIGEN3_INCLUDE_DIR "$ENV{EIGEN_DIR}" )
if (NOT EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 3.3 REQUIRED)
endif()
include_directories( ${EIGEN3_INCLUDE_DIR} )

# compiler flags
include(cmake/flags.cmake)

# utilities for cmake
include(cmake/utilities.cmake)

# set the CPU and GPU threading options
include(cmake/threading.cmake)

set( avro_EXTERNAL_LIBRARIES )

# optional libraries
list( APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake )
if (AVRO_WITH_MPI)
  get_dependency(MPI)
  set( MPI_INCLUDE_DIRS ${MPI_C_INCLUDE_PATH} )
  include_directories( SYSTEM ${MPI_INCLUDE_DIRS} )
  foreach (MPI_LIB ${MPI_LIBRARIES})
    get_filename_component( MPI_LIBRARY_DIR ${MPI_LIB} PATH )
    set_property( DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY LINK_DIRECTORIES ${MPI_LIBRARY_DIR} )
  endforeach()
  list( APPEND avro_EXTERNAL_LIBRARIES ${MPI_LIBRARIES} )

  get_dependency(PARMETIS)
  include_directories( SYSTEM ${PARMETIS_INCLUDE_DIRS} )
  list( APPEND avro_EXTERNAL_LIBRARIES ${PARMETIS_LIBRARIES} )
else()
  macro( RUN_WITH_MPIEXEC )
  endmacro()
endif()

# debug and coverage utilities
include(cmake/debug.cmake)
include(cmake/coverage.cmake)

set(AVRO_VERSION_MAJOR 2)
set(AVRO_VERSION_MINOR 1)
set(AVRO_VERSION_PATCH 0)
set(AVRO_VERSION ${AVRO_VERSION_MAJOR}.${AVRO_VERSION_MINOR}.${AVRO_VERSION_PATCH})

# Determine the current build date
string(TIMESTAMP avro_BUILD_DATE "%Y-%m-%d %H:%M:%S")
string(TIMESTAMP YEAR "%Y")

# Determine the current build number
# This is set by Jenkins in environment variable BUILD_NUMBER
set(avro_BUILD_NUMBER $ENV{BUILD_NUMBER})

# RPATH (where executables find the .so / DLLs)
# - Enables RPATH support for MACOSX
# - Makes RPATH of dynamic libraries and executable point to the directory
#   where libraries are installed.
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

add_definitions( -DGLIBCXX_USE_CXX11_ABI=0 )
#add_definitions( -DAVRO_SOURCE_DIR="${CMAKE_SOURCE_DIR}" )

include_directories( ${CMAKE_BINARY_DIR}/include )
add_subdirectory(src/third_party)
add_subdirectory(src)
add_subdirectory(test)

# this needs to happen here so the "avro" target is correct
add_subdirectory(src/bin)
add_subdirectory(src/examples)

function(config_bool cmake_var header_var )
  if (${cmake_var})
    set( ${header_var} 1 PARENT_SCOPE )
  else()
    set( ${header_var} 0 PARENT_SCOPE )
  endif()
endfunction()

# generate a configuration file
config_bool( AVRO_WITH_MPI CMAKE_AVRO_WITH_MPI )
config_bool( AVRO_HEADLESS_GRAPHICS CMAKE_AVRO_HEADLESS_GRAPHICS )
config_bool( AVRO_NO_ESP CMAKE_AVRO_NO_ESP )
config_bool( AVRO_WITH_GL CMAKE_AVRO_WITH_GL )
set( CMAKE_AVRO_CPU_THREAD_MANAGER ${AVRO_CPU_THREAD_MANAGER} )
configure_file( "${CMAKE_SOURCE_DIR}/src/include/avro_config.h.in" "${CMAKE_BINARY_DIR}/include/avro_config.h" )
