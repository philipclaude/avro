if (URSA_WITH_GL)

  include(FindOpenGL CONFIG)


  find_package(glfw QUIET)

  if (glf3_FOUND)
    message( "found glfw3 on the system: ${GLFW3_LIB}")
    set( URSA_BUILTIN_GLFW3 false )
  else()
    message(" compiling built-in glfw3")
    set( URSA_BUILTIN_GLFW3 true )
  endif()

else()
  message( "defaulting to only supporting webviewer")
endif()
