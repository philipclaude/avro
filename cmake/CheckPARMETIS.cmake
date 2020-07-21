include(CheckCXXSourceCompiles)

# make sure that the PARMETIS library is available
SET(CMAKE_REQUIRED_LIBRARIES ${PARMETIS_LIBRARIES} ${MPI_LIBRARIES})
SET(CMAKE_REQUIRED_INCLUDES ${PARMETIS_INCLUDE_DIRS} ${MPI_INCLUDE_DIRS})
IF( CMAKE_BUILD_TYPE )
  STRING( TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE )
  SET(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TYPE}} ${CMAKE_EXE_LINKER_FLAGS_${BUILD_TYPE}}")
ELSE()
  SET(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX_FLAGS}")
ENDIF()

IF( DEFINED PARMETIS_TEST AND (NOT PARMETIS_TEST OR NOT (PARMETIS_TEST_LIBRARIES STREQUAL PARMETIS_LIBRARIES AND PARMETIS_INCLUDE_DIRS STREQUAL PARMETIS_TEST_INCLUDE_DIRS)) )
  UNSET( PARMETIS_TEST CACHE )
  UNSET( PARMETIS_TEST_LIBRARIES CACHE )
  UNSET( PARMETIS_TEST_INCLUDE_DIRS CACHE )
ENDIF()

CHECK_CXX_SOURCE_COMPILES(
"
#include <parmetis.h>

#if PARMETIS_MAJOR_VERSION == 3 && PARMETIS_MINOR_VERSION < 2
#error \"Parmetis versions less than 3.2 have a memory allocation bug\"
#endif

#if PARMETIS_MAJOR_VERSION == 3
#define PARM_INT  idxtype
#define PARM_REAL float
#elif  PARMETIS_MAJOR_VERSION == 4
#define PARM_INT  idx_t
#define PARM_REAL real_t
#else
#error \"Unknown version of parmetis\"
#endif

int main(int argc, char *argv[])
{
  PARM_INT nparts, options[10];
  PARM_INT *part=NULL;
  PARM_REAL *tpwgts=NULL, *ubvec=NULL;
  MPI_Comm comm;
  PARM_INT numflag=0, wgtflag=0, edgecut=0;

  PARM_INT ncon=0;
  PARM_INT *xadj=NULL;          /* Pointers to the locally stored vertices */
  PARM_INT *vwgt=NULL;          /* Vertex weights */
  PARM_INT *adjncy=NULL;        /* Array that stores the adjacency lists of nvtxs */
  PARM_INT *adjwgt=NULL;        /* Array that stores the weights of the adjacency lists */
  PARM_INT *vtxdist=NULL;       /* Distribution of vertices */

  options[0] = 1; /* Turn on timing information */

  MPI_Init(&argc, &argv);
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  wgtflag = 3;
  ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt,
                       adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec,
                       options, &edgecut, part, &comm);

  MPI_Finalize();
}

" PARMETIS_TEST)

IF(PARMETIS_TEST)
  SET(PARMETIS_TEST_FAIL FALSE)

  # Save off the current state in case the library or inlcude is changed
  SET( PARMETIS_TEST_LIBRARIES ${PARMETIS_LIBRARIES} CACHE INTERNAL "PARMETIS libraries used for testing" FORCE)
  SET( PARMETIS_TEST_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIRS} CACHE INTERNAL "PARMETIS inlude used for testing" FORCE)
ELSE()

  MESSAGE(STATUS "" )
  MESSAGE(STATUS "====================================================================" )
  MESSAGE(STATUS " PARMETIS library test failed. See CMakeFiles/CMakeError.log for more details." )
  MESSAGE(STATUS "" )
  MESSAGE(STATUS " PARMETIS_LIBRARIES:" )
  FOREACH(LIBNAME ${PARMETIS_LIBRARIES})
  MESSAGE(STATUS "                    ${LIBNAME}" )
  ENDFOREACH()
  MESSAGE(STATUS "" )
  MESSAGE(STATUS " PARMETIS_INCLUDE_DIRS:" )
  FOREACH(INCNAME ${PARMETIS_INCLUDE_DIRS})
  MESSAGE(STATUS "                       ${INCNAME}" )
  ENDFOREACH()
  MESSAGE(STATUS "" )
  MESSAGE(STATUS " SANS will attempt to compile parmetis for you to resolve this." )
  MESSAGE(STATUS "====================================================================" )
  MESSAGE(STATUS "" )

  SET(PARMETIS_TEST_FAIL TRUE)
ENDIF()
