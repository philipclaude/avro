if (URSA_WITH_FPG)
	add_subdirectory(PCK/fpg)
endif()

#add_source_directories( TET_SOURCES "tetgen source" tetgen1.5.0 )
#set_source_files_properties( ${TET_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
#set_source_files_properties( ${TET_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY )

#add_definitions( -DTRILIBRARY )
#add_source_directories( TRIANGLE_SOURCES "triangle source" triangle )
#set_source_files_properties( ${TRIANGLE_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )

add_source_directories( TE_SOURCES "tiny expr source" tinyexpr )
set_source_files_properties( ${TE_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
include_directories( tinyexpr )

add_source_directories( JSON_SOURCES "json source" json )
set_source_files_properties( ${JSON_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
include_directories( json )

set( GLFW_SOURCES )
if (URSA_BUILTIN_GLFW)
	add_subdirectory( glfw )
endif()

set( GLAD_SOURCES )
if (URSA_WITH_GL)
	add_source_directories( GLAD_SOURCES "glad source" glad/src )
	include_directories( glad/include )
endif()

#add_source_directories( GMF_SOURCES "gmf source" libMeshb/libMeshb7 )
#set_source_files_properties( ${GMF_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" TRANSMESH )
#set_source_files_properties( ${GMF_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
#include_directories( libMeshb/libMeshb7 )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/tools )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/dynamic )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/dynamic/Eigen )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/dynamic/lapack )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/dynamic/MatMul )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/static )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/static/Eigen )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/static/MatMul )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/dense/static/SVD )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/block )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" numpack/block/tools )

add_library(ursa_third_party OBJECT ${SOURCES} ${TET_SOURCES} ${GMF_SOURCES} ${TE_SOURCES} ${TRIANGLE_SOURCES} ${NUMPACK_SOURCES} ${JSON_SOURCES} ${GLFW_SOURCES} ${GLAD_SOURCES} )
