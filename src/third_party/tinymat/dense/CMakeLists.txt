
INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( dense_SRC
     dynamic/MatMul/MatMul_Native_Real.cpp
     dynamic/MatMul/MatMul_Native_MatrixD_Instantiate.cpp
     dynamic/MatMul/MatMul_Native_MatrixS_Instantiate.cpp
     dynamic/MatMul/MatMul_BLAS_float.cpp
     dynamic/MatMul/MatMul_BLAS_double.cpp
     dynamic/MatrixD_Decompose_LU_Real.cpp
     dynamic/MatrixD_Decompose_LU_MatrixD_Instantiate.cpp
     dynamic/MatrixD_Decompose_LU_MatrixS_Instantiate.cpp
     dynamic/MatrixD_Decompose_LUP_Real.cpp
     dynamic/MatrixD_Decompose_LUP_MatrixD_Instantiate.cpp
     dynamic/MatrixD_Decompose_LUP_MatrixS_Instantiate.cpp
     dynamic/MatrixD_Decompose_QR_Real.cpp
     dynamic/MatrixD_Decompose_QR_MatrixS_Instantiate.cpp
     dynamic/MatrixD_InverseLU_Real_Instantiate.cpp
     dynamic/MatrixD_InverseLU_MatrixD_Instantiate.cpp
     dynamic/MatrixD_InverseLU_MatrixS_Instantiate.cpp
     dynamic/MatrixD_InverseLUP_Real_Instantiate.cpp
     dynamic/MatrixD_InverseLUP_MatrixD_Instantiate.cpp
     dynamic/MatrixD_InverseLUP_MatrixS_Instantiate.cpp
     dynamic/MatrixD_InverseQR_Real.cpp
     dynamic/MatrixD_InverseQR_MatrixS_Instantiate.cpp
     dynamic/MatrixD_Mul.cpp
     dynamic/MatrixD.cpp

     dynamic/ElementaryReflector_Real.cpp
     dynamic/ElementaryReflector_MatrixS_Instantiate.cpp

     dynamic/lapack/Eigen_double.cpp
     dynamic/lapack/Eigen_float.cpp

     static/MatMul/MatrixS_MatMul_Native_Int_Instantiate.cpp
     static/MatMul/MatrixS_MatMul_Native_Real_Instantiate.cpp
     static/MatMul/MatrixS_MatMul_Native_MatrixS_Instantiate.cpp
     static/MatMul/MatrixS_MatMul_Native_Surreal_Instantiate.cpp
     static/MatrixS_Decompose_LU_Real_Instantiate.cpp
     static/MatrixS_Decompose_LU_Surreal_Instantiate.cpp
     static/MatrixS_Decompose_LUP_Real_Instantiate.cpp
     static/MatrixS_Decompose_LUP_Surreal_Instantiate.cpp
     static/MatrixS_Decompose_QR_Real_Instantiate.cpp
     static/MatrixS_InverseLU_Real_Instantiate.cpp
     static/MatrixS_InverseLU_Surreal_Instantiate.cpp
     static/MatrixS_InverseLUP_Real_Instantiate.cpp
     static/MatrixS_InverseLUP_Surreal_Instantiate.cpp
     static/MatrixS_InverseQR_Real_Instantiate.cpp

     static/MatrixSymS_Decompose_Cholesky_Real_Instantiate.cpp
     static/MatrixSymS_Decompose_LDLT_Real_Instantiate.cpp
     static/MatrixSymS_InverseCholesky_Real_Instantiate.cpp
     static/MatrixSymS_InverseLDLT_Real_Instantiate.cpp

     static/ElementaryReflector_Real_Instantiate.cpp

     static/Eigen/MatrixSymS_Eigen_1x1.cpp
     static/Eigen/MatrixSymS_Eigen_2x2.cpp
     static/Eigen/MatrixSymS_Eigen_3x3.cpp
     static/Eigen/MatrixSymS_Eigen_4x4.cpp

     static/SVD/MatrixS_SVD_1x1.cpp
     static/SVD/MatrixS_SVD_2x2.cpp
     static/SVD/MatrixS_SVD_3x3.cpp
     static/SVD/MatrixS_SVD_4x4.cpp
     static/SVD/MatrixS_SVD_lapack_double.cpp
     static/SVD/MatrixS_SVD_lapack_float.cpp

     tools/Matrix_Util.cpp
     tools/SingularException.cpp
  )

ADD_LIBRARY( denseLib STATIC ${dense_SRC} )

#Create the vera targest for this library
ADD_VERA_CHECKS_RECURSE( denseLib *.h *.cpp )
ADD_HEADER_COMPILE_CHECK_RECURSE( denseLib *.h )
ADD_CPPCHECK( denseLib ${dense_SRC} )

IF( DEFINED LAPACK_DEPENDS )
  ADD_DEPENDENCIES( denseLib ${LAPACK_DEPENDS} )
  IF( USE_HEADERCHECK )
    ADD_DEPENDENCIES( denseLib_headercheck ${LAPACK_DEPENDS} )
  ENDIF()
ENDIF()

IF( DEFINED CBLAS_DEPENDS )
  ADD_DEPENDENCIES( denseLib ${CBLAS_DEPENDS} )
  IF( USE_HEADERCHECK )
    ADD_DEPENDENCIES( denseLib_headercheck ${CBLAS_DEPENDS} )
  ENDIF()
ENDIF()

IF( CMAKE_COMPILER_IS_GNUCXX )
  # Seems to be a bug in gcc array bounds for this file...
  SET_SOURCE_FILES_PROPERTIES( static/MatrixS_InverseLUP_Real_Instantiate.cpp PROPERTIES COMPILE_FLAGS "-Wno-array-bounds" )
ENDIF()


# USED FOR DEBUGGING
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()
