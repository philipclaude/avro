
INCLUDE( ${CMAKE_SOURCE_DIR}/CMakeInclude/ForceOutOfSource.cmake ) #This must be the first thing included

SET( SPARSELINALG_SRC
    SparseLinAlg_LinearSolver.cpp
    
    ScalarMatrix_CRS_int.cpp
    
    WriteMatrixMarketFile.cpp
    WritePlainVector.cpp
    
    Direct/UMFPACK/UMFPACKSolver.cpp
    Direct/UMFPACK/UMFPACKSolver_CRS_double.cpp
    Direct/UMFPACK/UMFPACKSolver_CRS_MatrixD.cpp
    Direct/UMFPACK/UMFPACKSolver_CRS_MatrixS_Instantiate.cpp
    Direct/UMFPACK/UMFPACKSolver_CRS_MatrixD_MatrixS_Instantiate.cpp
    Direct/UMFPACK/UMFPACKSolver_MatrixD_CRS_Real.cpp
    Direct/UMFPACK/UMFPACKSolver_MatrixD_CRS_MatrixD.cpp
    Direct/UMFPACK/UMFPACKSolver_MatrixD_CRS_MatrixS_Instantiate.cpp
    Direct/UMFPACK/UMFPACKSolver_MatrixD_CRS_MatrixD_MatrixS_Instantiate.cpp
    Direct/UMFPACK/UMFPACKSolver_MatrixBlock_2x2_Instantiate.cpp
    Direct/UMFPACK/UMFPACKSolver_MatrixBlock_3x3_Instantiate.cpp
    Direct/UMFPACK/UMFPACKSolver_MatrixBlock_4x4_Instantiate.cpp
    Direct/UMFPACK/ScalarMatrix_CRS_UF_INT.cpp

    Krylov/FGMRES.cpp
    Krylov/FGMRES_CRS_Real.cpp
    Krylov/FGMRES_CRS_MatrixD.cpp
    Krylov/FGMRES_CRS_MatrixS_Instantiate.cpp
    Krylov/FGMRES_CRS_MatrixD_MatrixS_Instantiate.cpp
    Krylov/FGMRES_MatrixD_CRS_Real.cpp
    Krylov/FGMRES_MatrixD_CRS_MatrixD.cpp
    Krylov/FGMRES_MatrixD_CRS_MatrixS_Instantiate.cpp
    Krylov/FGMRES_MatrixD_CRS_MatrixD_MatrixS_Instantiate.cpp
    Krylov/FGMRES_MatrixBlock_2x2_Instantiate.cpp
    Krylov/FGMRES_MatrixBlock_4x4_Instantiate.cpp
    
    Preconditioners/Identity.cpp

    Preconditioners/LU_SGS.cpp
    Preconditioners/LU_SGS_Real_Instantiate.cpp
    Preconditioners/LU_SGS_MatrixD_Instantiate.cpp
    Preconditioners/LU_SGS_MatrixS_Instantiate.cpp
    Preconditioners/LU_SGS_MatrixD_MatrixS_Instantiate.cpp
    Preconditioners/LU_SGS_MatrixD_Sparse_MatrixS_Instantiate.cpp
    Preconditioners/LU_SGS_MatrixBlock_2x2_Instantiate.cpp

    Preconditioners/ILU0.cpp
    Preconditioners/ILU0_Real_Instantiate.cpp
    
    tools/UpperTriangMatrix.cpp
  )
 
# MKL is optional
IF( USE_MKL )

LIST( APPEND SPARSELINALG_SRC

    Direct/MKL_PARDISO/MKL_PARDISOSolver.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_CRS_double.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_CRS_MatrixD.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_CRS_MatrixS_Instantiate.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_CRS_MatrixD_MatrixS_Instantiate.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_MatrixD_CRS_Real.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_MatrixD_CRS_MatrixD.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_MatrixD_CRS_MatrixS_Instantiate.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_MatrixD_CRS_MatrixD_MatrixS_Instantiate.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_MatrixBlock_2x2_Instantiate.cpp
    Direct/MKL_PARDISO/MKL_PARDISOSolver_MatrixBlock_4x4_Instantiate.cpp
    Direct/MKL_PARDISO/ScalarMatrix_CRS_MKL_INT.cpp

  )

ENDIF()

# PETSc is optional
IF( USE_PETSC )

LIST( APPEND SPARSELINALG_SRC

    PETSc/PETScSolver.cpp
    PETSc/PETScSolver_CRS_double.cpp
    PETSc/PETScSolver_CRS_MatrixD.cpp
    PETSc/PETScSolver_CRS_MatrixS_Instantiate.cpp
    PETSc/PETScSolver_CRS_MatrixD_MatrixS_Instantiate.cpp
    PETSc/PETScSolver_MatrixD_CRS_Real_Instantiate.cpp
    PETSc/PETScSolver_MatrixD_CRS_MatrixD.cpp
    PETSc/PETScSolver_MatrixD_CRS_MatrixS_Instantiate.cpp
    PETSc/PETScSolver_MatrixD_CRS_MatrixD_MatrixS_Instantiate.cpp
    PETSc/PETScSolver_MatrixBlock_2x2_Instantiate.cpp
    PETSc/PETScSolver_MatrixBlock_3x3_Instantiate.cpp
    PETSc/PETScSolver_MatrixBlock_4x4_Instantiate.cpp
    PETSc/MDFOrdering.cpp
  )

ENDIF()

ADD_LIBRARY( SparseLinAlgLib STATIC ${SPARSELINALG_SRC} )

# Need a list of files for header check and vera
SET( SPARSELINALG_CHECK_H *.h Direct/UMFPACK/*.h
                              Krylov/*.h
                              Preconditioners/*.h )

SET( SPARSELINALG_CHECK_CPP *.cpp Direct/UMFPACK/*.cpp
                                  Krylov/*.cpp
                                  Preconditioners/*cpp )

# MKL is optional
IF( USE_MKL )
  LIST( APPEND SPARSELINALG_CHECK_H   Direct/MKL_PARDISO/*.h   )
  LIST( APPEND SPARSELINALG_CHECK_CPP Direct/MKL_PARDISO/*.cpp )
ENDIF()

# PETSc is optional
IF( USE_PETSC )
  LIST( APPEND SPARSELINALG_CHECK_H   PETSc/*.h   )
  LIST( APPEND SPARSELINALG_CHECK_CPP PETSc/*.cpp )
ENDIF()

#Create the vera targest for this library
ADD_VERA_CHECKS( SparseLinAlgLib ${SPARSELINALG_CHECK_H} ${SPARSELINALG_CHECK_CPP} )
ADD_HEADER_COMPILE_CHECK( SparseLinAlgLib ${SPARSELINALG_CHECK_H} )
ADD_CPPCHECK( SparseLinAlgLib ${SPARSELINALG_SRC} )

IF( DEFINED PETSC_DEPENDS )
  ADD_DEPENDENCIES( SparseLinAlgLib ${PETSC_DEPENDS} )
  IF( USE_HEADERCHECK )
    ADD_DEPENDENCIES( SparseLinAlgLib_headercheck ${PETSC_DEPENDS} )
  ENDIF()
ENDIF()

IF( DEFINED SUITESPARSE_DEPENDS )
  ADD_DEPENDENCIES( SparseLinAlgLib ${SUITESPARSE_DEPENDS} )
  IF( USE_HEADERCHECK )
    ADD_DEPENDENCIES( SparseLinAlgLib_headercheck ${SUITESPARSE_DEPENDS} )
  ENDIF()
ENDIF()
