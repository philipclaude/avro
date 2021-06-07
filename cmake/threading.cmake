# always look for OpenMP
if (avro_CPU_THREAD_MANAGER STREQUAL "openmp")
  include(cmake/FindOpenMP.cmake)
  if (NOT OPENMP_FOUND)

    # cannot use OpenMP as thread manager if it isn't found!
    set( avro_WITH_OPENMP false )
    set( avro_CPU_THREAD_MANAGER "cpp11" )
    message( STATUS "OpenMP not found: defaulting to c++11 threads.")

  else()

    set( avro_CPU_THREAD_MANAGER "openmp" )
    set( avro_WITH_OPENMP true )
    message( STATUS "CPU thread manager = OpenMP.")
    add_definitions(-Davro_CPU_THREAD_MANAGER_OPENMP)
    
  endif()

else()
  set( avro_WITH_OPENMP false )
endif()

if (avro_CPU_THREAD_MANAGER STREQUAL "cpp11")
  # set the CPU thread manager as c++11 threads
  message( STATUS "CPU thread manager = c++11.")
  add_definitions(-Davro_CPU_THREAD_MANAGER_CPP)
endif()

if (avro_CPU_THREAD_MANAGER STREQUAL "pthread")
  # set the CPU thread manager as pthreads
  message( STATUS "CPU thread manager = pthread.")
  add_definitions(-Davro_CPU_THREAD_MANAGER_PTHREAD)
endif()

if (avro_CPU_THREAD_MANAGER STREQUAL "emp")
  # set the CPU thread manager as EMP threads
  message( STATUS "CPU thread manager = EMP (part of EngSketchPad).")
  add_definitions(-Davro_CPU_THREAD_MANAGER_EMP)
endif()

# determine the GPU threading environment
if (avro_GPU_THREAD_MANAGER STREQUAL "cuda")

  # check if we have cuda
  include(cmake/FindCUDA.cmake)
  if (NOT CUDA_FOUND)
    message(FATAL_ERROR "CUDA was not found")
  endif()
  add_definitions(-Davro_GPU_THREAD_MANAGER_CUDA)
  set(avro_WITH_CUDA true )
elseif( avro_GPU_THREADING STREQUAL "opencl" )

  # check if we have cuda
  include(cmake/FindOpenCLcmake)
  if (NOT OPENCL_FOUND)
    message(FATAL_ERROR "OpenCL was not found")
  endif()
  add_definitions(-Davro_GPU_THREAD_MANAGER_OPENCL)
  set(avro_WITH_OPENCL true)

else()
  message( STATUS "GPU threading not enabled. GPU thread manager = serial." )
endif()
