INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( block_SRC

	MatrixBlock_2x2_Decompose_LU.cpp
	MatrixBlock_2x2_Decompose_LUP.cpp
	MatrixBlock_2x2_InverseLU.cpp
	MatrixBlock_2x2_InverseLUP.cpp
	
   )

ADD_LIBRARY( blockLib STATIC ${block_SRC} )

#Create the vera targets for this library
ADD_VERA_CHECKS( blockLib *.h *.cpp )
ADD_HEADER_COMPILE_CHECK( blockLib *.h )
ADD_CPPCHECK( blockLib ${block_SRC} )