INCLUDE(ExternalProject)

SET( PARMETIS_VERSION 4.0.3 )

SET( PARMETIS_URL ${CMAKE_SOURCE_DIR}/external/parmetis-${PARMETIS_VERSION}.tar.gz )
MESSAGE(STATUS "Checking for ${PARMETIS_URL}")
IF( NOT EXISTS ${PARMETIS_URL} )
  SET( PARMETIS_URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-${PARMETIS_VERSION}.tar.gz )
  MESSAGE(STATUS "Downloading and compiling ParMETIS from ${PARMETIS_URL}")
ELSE()
  MESSAGE(STATUS "Using ${PARMETIS_URL}")
ENDIF()

SET(PARMETIS_INSTALL_DIR ${CMAKE_BINARY_DIR}/libparmetis-prefix/install)
SET(PARMETIS_SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/parmetis-${PARMETIS_VERSION} )

#SET( PARMETIS_SOURCE_DIR ${CMAKE_BINARY_DIR}/libparmetis-prefix/ParMETIS-${PARMETIS_VERSION} )
#SET( PARMETIS_INSTALL_DIR ${PARMETIS_SOURCE_DIR}/install/ )

# Patches from PETSc
# https://bitbucket.org/petsc/pkg-parmetis/commits/82409d68aa1d6cbc70740d0f35024aae17f7d5cb
# https://bitbucket.org/petsc/pkg-parmetis/commits/1c1a9fd0f408dc4d42c57f5c3ee6ace411eb222b/
# http://glaros.dtc.umn.edu/gkhome/node/837

FIND_PROGRAM( PATCH patch REQUIRED )

STRING (REPLACE ";" " " PARMETIS_LIBS "${MPI_C_LIBRARIES}")
STRING (REPLACE ";" " -I" PARMETIS_CFLAGS "-I${MPI_C_INCLUDE_PATH} -w")

STRING( TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE )
IF( BUILD_TYPE MATCHES "DEBUG" )
  # NOTE: -DGDB=1 causes the "warning: label ‘CleanUpAndExit’ defined but not used [-Wunused-label]" to be flagged as an error (-Werror)
  SET( PARMETIS_CMAKE_ARGS  -DASSERT=1 -DASSERT2=1 -DDEBUG=1)
ENDIF()

IF( APPLE )
  SET( PARMETS_SHARED_LINKER_FLAGS "-undefined dynamic_lookup" )
ENDIF()

SET( PARMETIS_CMAKE_ARGS ${PARMETIS_CMAKE_ARGS}
                         -DSHARED=1
                         -DCMAKE_INSTALL_PREFIX=${PARMETIS_INSTALL_DIR}
                         -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                         -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                         -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                         -DCMAKE_C_FLAGS=${PARMETIS_CFLAGS}
                         -DCMAKE_SHARED_LINKER_FLAGS=${PARMETS_SHARED_LINKER_FLAGS}
                         -DMETIS_INSTALL=ON
   )

SET( PARMETIS_CMAKE_CACHE_ARGS -DGKLIB_PATH:STRING=${PARMETIS_SOURCE_DIR}/metis/GKlib
                               -DMETIS_PATH:STRING=${PARMETIS_SOURCE_DIR}/metis
                               -DCMAKE_C_FLAGS_MEMCHECK:STRING=${CMAKE_C_FLAGS_MEMCHECK}
                               -DCMAKE_SHARED_LINKER_FLAGS_MEMCHECK:STRING=${CMAKE_SHARED_LINKER_FLAGS_MEMCHECK}
   )


ExternalProject_Add(
  libparmetis
  DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/external/
  URL_MD5 f69c479586bf6bb7aff6a9bc0c739628
  URL ${PARMETIS_URL}
  SOURCE_DIR ${PARMETIS_SOURCE_DIR}
  PATCH_COMMAND ${PATCH} -p1 -i ${CMAKE_SOURCE_DIR}/cmake/patch/parmetis-xyzpart.patch
  COMMAND       ${PATCH} -p0 -i ${CMAKE_SOURCE_DIR}/cmake/patch/parmetis-CMakeLists.patch
  COMMAND       ${PATCH} -p0 -i ${CMAKE_SOURCE_DIR}/cmake/patch/parmetis-metis-CMakeLists.patch
  COMMAND       ${PATCH} -p0 -i ${CMAKE_SOURCE_DIR}/cmake/patch/parmetis-libparmetis-defs.patch
  CMAKE_CACHE_ARGS ${PARMETIS_CMAKE_CACHE_ARGS}
  CMAKE_ARGS ${PARMETIS_CMAKE_ARGS}
  )
#  CONFIGURE_COMMAND $(MAKE) config prefix=${PARMETIS_INSTALL_DIR} cc=${MPI_C_COMPILER} cxx=${MPI_CXX_COMPILER} ${PARMETIS_ARGS}
#  BUILD_COMMAND $(MAKE) MAKEFLAGS=

SET_PROPERTY( TARGET libparmetis PROPERTY FOLDER "Externals")

SET( PARMEITS_ARCH ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR} )

SET(PARMETIS_INCLUDE_DIRS ${PARMETIS_INSTALL_DIR}/include)
SET(PARMETIS_LIBRARIES ${PARMETIS_INSTALL_DIR}/lib/libparmetis${CMAKE_SHARED_LIBRARY_SUFFIX}
                       ${PARMETIS_INSTALL_DIR}/lib/libmetis${CMAKE_SHARED_LIBRARY_SUFFIX}  )
SET(PARMETIS_DEPENDS libparmetis)
SET(PARMETIS_USE_INTERNAL TRUE)

SET_PROPERTY( DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY
              LINK_DIRECTORIES ${PARMETIS_INSTALL_DIR}/lib )

LIST(APPEND SANS_EXTERNAL_DEPENDS ${PARMETIS_DEPENDS})
