if (URSA_WITH_GL)

  include(FindOpenGL CONFIG)


  find_package(glfw QUIET)

  if (glf3_FOUND)
    message( STATUS "found glfw3 on the system: ${GLFW3_LIB}")
    set( URSA_BUILTIN_GLFW false )
  else()
    message(STATUS "compiling built-in glfw3")
    set( URSA_BUILTIN_GLFW true )
  endif()

  if (APPLE)
    find_library(COCOA_LIBRARY Cocoa)
  endif()
else()
  message( "defaulting to only supporting webviewer")
endif()
