cmake_minimum_required(VERSION 2.8.8)

# this allows setting the version number
if (POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif()

# setup libraries and build type
include(cmake/setup.cmake)
include(cmake/compiler.cmake)
project(luma VERSION 1.0)



if (POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)
endif()

# default options
option(luma_STATIC_LINK_CXX "statically link to c++ libraries (probably because you're using voyager)" false)
option(luma_SMART_PTR "use internal smart pointer (smartypants) for memory management" false )
option(luma_WITH_MPI "compile for distributed computing support" false)
option(luma_RUN_WITH_MPI "run with mpi. might be difficult to debug with mpirun" false)
option(luma_WITH_FPG "compile FPG needed to make new predicates" false)
set(luma_TEST_DATA "data/" CACHE STRING "location for luma test files relative to test/" )
set(luma_WITH_GL "option to compile OpenGL-based visualizer (along with WebGL-based one)" true)


# PkgConfig is used by many Finds/*.cmake
FIND_PACKAGE(PkgConfig)

set(luma_CPU_THREAD_MANAGER "cpp11" CACHE STRING "type of cpu threading environment to use [openmp, cpp11, pthread, emp]")
set(luma_GPU_THREAD_MANAGER "none" CACHE STRING "type of gpu threading environment to use [opencl,cuda]")

# schedule dynamic doesn't seem to be supported by intel
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  set(luma_CPU_THREAD_MANAGER "cpp11")
endif()

# luma requires the c++11 standard
set(CMAKE_CXX_STANDARD 11)

# option to use internal smart pointer implementation
if (luma_SMART_PTR)
	add_definitions(-Dluma_SMART_PTR)
endif()

# required libraries
include(cmake/FindLAPACK.cmake)
include(cmake/FindESP.cmake)
include(cmake/FindNLOPT.cmake)

# optional libraries
if (luma_WITH_MPI)
  include(cmake/FindMPI.cmake)
  include(cmake/FindPARMETIS.cmake)
endif()

include(cmake/FindGraphics.cmake)

# compiler flags
include(cmake/flags.cmake)

# utilities for cmake
include(cmake/utilities.cmake)

# set the CPU and GPU threading options
include(cmake/threading.cmake)

if (NOT luma_RUN_WITH_MPI)
    add_definitions( -DNO_MPIRUN )
endif()

set( luma_EXTERNAL_LIBRARIES )

# debug and coverage utilities
include(cmake/debug.cmake)
include(cmake/coverage.cmake)

set(luma_VERSION_MAJOR 1)
set(luma_VERSION_MINOR 0)
set(luma_VERSION_PATCH 0)
set(luma_VERSION ${luma_VERSION_MAJOR}.${luma_VERSION_MINOR}.${luma_VERSION_PATCH})

# Determine the current build date
string(TIMESTAMP luma_BUILD_DATE "%Y-%m-%d %H:%M:%S")
string(TIMESTAMP YEAR "%Y")

# Determine the current build number
# This is set by Jenkins in environment variable BUILD_NUMBER
set(luma_BUILD_NUMBER $ENV{BUILD_NUMBER})

# RPATH (where executables find the .so / DLLs)
# - Enables RPATH support for MACOSX
# - Makes RPATH of dynamic libraries and executable point to the directory
#   where libraries are installed.
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

add_definitions( -DGLIBCXX_USE_CXX11_ABI=0 )

add_subdirectory(src/third_party)
add_subdirectory(src)
add_subdirectory(test)

# this needs to happen here so the "luma" target is correct
add_subdirectory(src/bin)
