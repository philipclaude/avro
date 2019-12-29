


if (NOT luma_USE_CPP11)
  if (luma_CPU_THREADING STREQUAL "cpp11")
    message( STATUS "cannot use cpp11 threads: setting emp" )
    set( luma_CPU_THREADING "emp" )
  endif()
  if (NOT luma_SMART_PTR)
    message( STATUS "using internal smart pointers since not cpp11 (shared_ptr)" )
    set(luma_SMART_PTR true)
  endif()
endif()

# always look for OpenMP
include(cmake/FindOpenMP.cmake)
if (NOT OPENMP_FOUND)

  message( STATUS "OpenMP not found" )
  set( luma_WITH_OPENMP FALSE )

  if (luma_CPU_THREAD_MANAGER STREQUAL "openmp")
    # cannot use OpenMP as thread manager if it isn't found!
    message( STATUS "OpenMP not found: defaulting to c++11 threads.")
    set( luma_CPU_THREAD_MANAGER "cpp11" )
  endif()

else()

  # OpenMP was found, so let's use it!
  add_definitions( -Dluma_WITH_OPENMP )
  set( luma_WITH_OPENMP TRUE )

  # if CPU thread manager is openmp
  if (luma_CPU_THREAD_MANAGER STREQUAL "openmp")
    message( STATUS "CPU thread manager = OpenMP.")
    add_definitions(-Dluma_CPU_THREAD_MANAGER_OPENMP)
  endif()

endif()

if (luma_CPU_THREAD_MANAGER STREQUAL "cpp11")
  # set the CPU thread manager as c++11 threads
  message( STATUS "CPU thread manager = c++11.")
  add_definitions(-Dluma_CPU_THREAD_MANAGER_CPP)
endif()

if (luma_CPU_THREAD_MANAGER STREQUAL "pthread")
  # set the CPU thread manager as pthreads
  message( STATUS "CPU thread manager = pthread.")
  add_definitions(-Dluma_CPU_THREAD_MANAGER_PTHREAD)
endif()

if (luma_CPU_THREADING STREQUAL "emp")
  # set the CPU thread manager as EMP threads
  message( STATUS "CPU thread manager = EMP (part of EngSketchPad).")
  add_definitions(-Dluma_CPU_THREAD_MANAGER_EMP)
endif()

# determine the GPU threading environment
if (luma_GPU_THREAD_MANAGER STREQUAL "cuda")

  # check if we have cuda
  include(cmake/FindCUDA.cmake)
  if (NOT CUDA_FOUND)
    message(FATAL_ERROR "CUDA was not found")
  endif()
  add_definitions(-Dluma_GPU_THREAD_MANAGER_CUDA)
  set(luma_WITH_CUDA true )
elseif( luma_GPU_THREADING STREQUAL "opencl" )

  # check if we have cuda
  include(cmake/FindOpenCLcmake)
  if (NOT OPENCL_FOUND)
    message(FATAL_ERROR "OpenCL was not found")
  endif()
  add_definitions(-Dluma_GPU_THREAD_MANAGER_OPENCL)
  set(luma_WITH_OPENCL true)

else()
  message( STATUS "GPU threading not enabled. GPU thread manager = serial." )
endif()
