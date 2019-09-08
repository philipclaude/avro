
INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( SURREAL_SRC
     SurrealD.cpp
  )

ADD_LIBRARY( SurrealLib STATIC ${SURREAL_SRC} )

#Create the vera targest for this library
ADD_VERA_CHECKS_RECURSE( SurrealLib *.h *.cpp )
ADD_HEADER_COMPILE_CHECK_RECURSE( SurrealLib *.h )
ADD_CPPCHECK( SurrealLib ${SURREAL_SRC} )