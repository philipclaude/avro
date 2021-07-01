if (avro_WITH_GL)

  include(FindOpenGL CONFIG)

  set(CMAKE_THREAD_LIBS_INIT "-lpthread")
  set(CMAKE_HAVE_THREADS_LIBRARY 1)
  set(CMAKE_USE_WIN32_THREADS_INIT 0)
  set(CMAKE_USE_PTHREADS_INIT 1)
  set(THREADS_PREFER_PTHREAD_FLAG ON)

  find_package(glfw QUIET)

  if (glfw3_FOUND)
    message( STATUS "found glfw3 on the system: ${GLFW3_LIB}")
    set( avro_BUILTIN_GLFW false )
  else()
    message(STATUS "compiling built-in glfw3")
    set( avro_BUILTIN_GLFW true )
  endif()

  if (APPLE)
    find_library(COCOA_LIBRARY Cocoa)
  endif()

  set( AVRO_WITH_GL true )
else()
  message( STATUS "only supporting WebViewer")
endif()
