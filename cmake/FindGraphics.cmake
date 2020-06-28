if (avro_WITH_GL)

  include(FindOpenGL CONFIG)


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
else()
  message( "defaulting to only supporting webviewer")
endif()

#find_package(GLM CONFIG REQUIRED)
#include_directories( ${GLM_INCLUDE_DIR} )
