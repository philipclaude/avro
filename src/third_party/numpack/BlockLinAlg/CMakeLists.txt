INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( BLOCKLINALG_SRC

	MatrixBlock_2x2_Decompose_LU.cpp
	MatrixBlock_2x2_Decompose_LUP.cpp
	MatrixBlock_2x2_InverseLU.cpp
	MatrixBlock_2x2_InverseLUP.cpp
	
   )

ADD_LIBRARY( BlockLinAlgLib STATIC ${BLOCKLINALG_SRC} )

#Create the vera targets for this library
ADD_VERA_CHECKS( BlockLinAlgLib *.h *.cpp )
ADD_HEADER_COMPILE_CHECK( BlockLinAlgLib *.h )
ADD_CPPCHECK( BlockLinAlgLib ${BLOCKLINALG_SRC} )