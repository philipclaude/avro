
INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( TOOLS_SRC
     SANSException.cpp
     BackTraceException.cpp
     demangle.cpp
     split_string.cpp
     system_call.cpp
  )

ADD_LIBRARY( toolsLib STATIC ${TOOLS_SRC} )
#ADD_DEPENDENCIES( tools libbacktrace ) libbacktrace leaks memory...

#Create the vera targest for this library
ADD_VERA_CHECKS_RECURSE( toolsLib *.h *.cpp )
ADD_HEADER_COMPILE_CHECK_RECURSE( toolsLib *.h )
ADD_CPPCHECK( toolsLib ${TOOLS_SRC} )

# USED FOR DEBUGGING
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()
