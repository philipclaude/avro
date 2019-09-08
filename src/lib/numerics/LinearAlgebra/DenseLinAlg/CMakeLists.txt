
INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( DENSELINALG_SRC
     DynamicSize/MatMul/MatMul_Native_Real.cpp
     DynamicSize/MatMul/MatMul_Native_MatrixD_Instantiate.cpp
     DynamicSize/MatMul/MatMul_Native_MatrixS_Instantiate.cpp
     DynamicSize/MatMul/MatMul_BLAS_float.cpp
     DynamicSize/MatMul/MatMul_BLAS_double.cpp
     DynamicSize/MatrixD_Decompose_LU_Real.cpp
     DynamicSize/MatrixD_Decompose_LU_MatrixD_Instantiate.cpp
     DynamicSize/MatrixD_Decompose_LU_MatrixS_Instantiate.cpp
     DynamicSize/MatrixD_Decompose_LUP_Real.cpp
     DynamicSize/MatrixD_Decompose_LUP_MatrixD_Instantiate.cpp
     DynamicSize/MatrixD_Decompose_LUP_MatrixS_Instantiate.cpp
     DynamicSize/MatrixD_Decompose_QR_Real.cpp
     DynamicSize/MatrixD_Decompose_QR_MatrixS_Instantiate.cpp
     DynamicSize/MatrixD_InverseLU_Real_Instantiate.cpp
     DynamicSize/MatrixD_InverseLU_MatrixD_Instantiate.cpp
     DynamicSize/MatrixD_InverseLU_MatrixS_Instantiate.cpp
     DynamicSize/MatrixD_InverseLUP_Real_Instantiate.cpp
     DynamicSize/MatrixD_InverseLUP_MatrixD_Instantiate.cpp
     DynamicSize/MatrixD_InverseLUP_MatrixS_Instantiate.cpp
     DynamicSize/MatrixD_InverseQR_Real.cpp
     DynamicSize/MatrixD_InverseQR_MatrixS_Instantiate.cpp
     DynamicSize/MatrixD_Mul.cpp
     DynamicSize/MatrixD.cpp

     DynamicSize/ElementaryReflector_Real.cpp
     DynamicSize/ElementaryReflector_MatrixS_Instantiate.cpp

     DynamicSize/lapack/Eigen_double.cpp
     DynamicSize/lapack/Eigen_float.cpp

     StaticSize/MatMul/MatrixS_MatMul_Native_Int_Instantiate.cpp
     StaticSize/MatMul/MatrixS_MatMul_Native_Real_Instantiate.cpp
     StaticSize/MatMul/MatrixS_MatMul_Native_MatrixS_Instantiate.cpp
     StaticSize/MatMul/MatrixS_MatMul_Native_Surreal_Instantiate.cpp
     StaticSize/MatrixS_Decompose_LU_Real_Instantiate.cpp
     StaticSize/MatrixS_Decompose_LU_Surreal_Instantiate.cpp
     StaticSize/MatrixS_Decompose_LUP_Real_Instantiate.cpp
     StaticSize/MatrixS_Decompose_LUP_Surreal_Instantiate.cpp
     StaticSize/MatrixS_Decompose_QR_Real_Instantiate.cpp
     StaticSize/MatrixS_InverseLU_Real_Instantiate.cpp
     StaticSize/MatrixS_InverseLU_Surreal_Instantiate.cpp
     StaticSize/MatrixS_InverseLUP_Real_Instantiate.cpp
     StaticSize/MatrixS_InverseLUP_Surreal_Instantiate.cpp
     StaticSize/MatrixS_InverseQR_Real_Instantiate.cpp

     StaticSize/MatrixSymS_Decompose_Cholesky_Real_Instantiate.cpp
     StaticSize/MatrixSymS_Decompose_LDLT_Real_Instantiate.cpp
     StaticSize/MatrixSymS_InverseCholesky_Real_Instantiate.cpp
     StaticSize/MatrixSymS_InverseLDLT_Real_Instantiate.cpp

     StaticSize/ElementaryReflector_Real_Instantiate.cpp

     StaticSize/Eigen/MatrixSymS_Eigen_1x1.cpp
     StaticSize/Eigen/MatrixSymS_Eigen_2x2.cpp
     StaticSize/Eigen/MatrixSymS_Eigen_3x3.cpp
     StaticSize/Eigen/MatrixSymS_Eigen_4x4.cpp

     StaticSize/SVD/MatrixS_SVD_1x1.cpp
     StaticSize/SVD/MatrixS_SVD_2x2.cpp
     StaticSize/SVD/MatrixS_SVD_3x3.cpp
     StaticSize/SVD/MatrixS_SVD_4x4.cpp
     StaticSize/SVD/MatrixS_SVD_lapack_double.cpp
     StaticSize/SVD/MatrixS_SVD_lapack_float.cpp

     tools/Matrix_Util.cpp
     tools/SingularException.cpp
  )

ADD_LIBRARY( DenseLinAlgLib STATIC ${DENSELINALG_SRC} )

#Create the vera targest for this library
ADD_VERA_CHECKS_RECURSE( DenseLinAlgLib *.h *.cpp )
ADD_HEADER_COMPILE_CHECK_RECURSE( DenseLinAlgLib *.h )
ADD_CPPCHECK( DenseLinAlgLib ${DENSELINALG_SRC} )

IF( DEFINED LAPACK_DEPENDS )
  ADD_DEPENDENCIES( DenseLinAlgLib ${LAPACK_DEPENDS} )
  IF( USE_HEADERCHECK )
    ADD_DEPENDENCIES( DenseLinAlgLib_headercheck ${LAPACK_DEPENDS} )
  ENDIF()
ENDIF()

IF( DEFINED CBLAS_DEPENDS )
  ADD_DEPENDENCIES( DenseLinAlgLib ${CBLAS_DEPENDS} )
  IF( USE_HEADERCHECK )
    ADD_DEPENDENCIES( DenseLinAlgLib_headercheck ${CBLAS_DEPENDS} )
  ENDIF()
ENDIF()

IF( CMAKE_COMPILER_IS_GNUCXX )
  # Seems to be a bug in gcc array bounds for this file...
  SET_SOURCE_FILES_PROPERTIES( StaticSize/MatrixS_InverseLUP_Real_Instantiate.cpp PROPERTIES COMPILE_FLAGS "-Wno-array-bounds" )
ENDIF()


# USED FOR DEBUGGING
#get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
#foreach(dir ${dirs})
#  message(STATUS "dir='${dir}'")
#endforeach()
