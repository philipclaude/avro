if (avro_WITH_FPG)
	add_subdirectory(PCK/fpg)
endif()

add_source_directories( TET_SOURCES "tetgen source" tetgen1.5.0 )
set_source_files_properties( ${TET_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
set_source_files_properties( ${TET_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY )

add_definitions( -DTRILIBRARY )
add_source_directories( TRIANGLE_SOURCES "triangle source" triangle )
set_source_files_properties( ${TRIANGLE_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )

add_source_directories( JSON_SOURCES "json source" json )
set_source_files_properties( ${JSON_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
include_directories( json )

if (avro_WITH_GL)
	if (avro_BUILTIN_GLFW)
		add_subdirectory( glfw )
		set( THIRD_PARTY_INCLUDE_DIRS ${THIRD_PARTY_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/third_party/glfw/include )
		include_directories( ${CMAKE_SOURCE_DIR}/src/third_party/glfw/include )
		set_target_properties( glfw PROPERTIES COMPILE_FLAGS "-w" )
	endif()

	add_source_directories( IMGUI_SOURCES "imgui source" imgui )
	add_source_directories( IMGUI_SOURCES "imgui source" imgui/GL )
	set_source_files_properties( ${IMGUI_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" IMGUI_IMPL_OPENGL_LOADER_GLAD )
	include_directories( imgui )


	set( GLAD_SOURCES )
	if (avro_WITH_GL)
		add_source_directories( GLAD_SOURCES "glad source" glad/src )
		set( THIRD_PARTY_INCLUDE_DIRS ${THIRD_PARTY_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/third_party/glad/include )
		include_directories( ${CMAKE_SOURCE_DIR}/src/third_party/glad/include )
	endif()
else()
	set( IMGUI_SOURCES )
endif()

# nearest neighbours
add_source_directories( NNSEARCH_SOURCES "nnsearch" nnsearch )

# libmeshb
add_source_directories( MESHB_SOURCES "gmf source" libmeshb )
set_source_files_properties( ${MESHB_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
include_directories( libmeshb )

# blossom
if (NOT avro_LITE)
	add_source_directories( BLOSSOM_SOURCES "blossom source" blossom5 )
	add_source_directories( BLOSSOM_SOURCES "blossom source" blossom5/GEOM )
	add_source_directories( BLOSSOM_SOURCES "blossom source" blossom5/MinCost )
	set_source_files_properties( ${BLOSSOM_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
	include_directories( blossom5 )
else()
	set( BLOSSOM_SOURCES )
endif()

# tinymat
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat )
add_source_directories(NUMPACK_SOURCES "Source Files\\tinymat" tinymat/tools )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/tools )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/dynamic )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/dynamic/Eigen )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/dynamic/lapack )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/dynamic/MatMul )

add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/static )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/static/Eigen )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/static/MatMul )
add_source_directories(NUMPACK_SOURCES "Source Files\\numerics" tinymat/dense/static/SVD )

include_directories(tinymat)

#add_source_directories( QHULL_SOURCES "QHULL" qhull/libqhull )
#include_directories( qhull )

add_source_directories(HLBFGS_SOURCES "HLBFGS" HLBFGS )
set_source_files_properties( ${HLBFGS_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )

add_source_directories(OPENNL_SOURCES "OpenNL" OpenNL )

add_source_directories(LIBLBFGS_SOURCES "LIBLBFGS" liblbfgs )
set_source_files_properties( ${LIBLBFGS_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )


# ESP
if (AVRO_NO_ESP)
	add_source_directories( WV_SOURCES "wv" wv )
	add_source_directories( EGLITE_SOURCES "egads" egads )
	set_source_files_properties( ${WV_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" LWS_NO_FORK )
	set_source_files_properties( ${WV_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
	set_source_files_properties( ${EGLITE_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" LITE )
	set_source_files_properties( ${EGLITE_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
	set( ESP_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/src/third_party/egads ${CMAKE_SOURCE_DIR}/src/third_party/wv )
	set( ESP_LIBRARIES )
endif()

set( THIRD_PARTY_INCLUDE_DIRS ${THIRD_PARTY_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS} )

add_library(avro_third_party OBJECT ${LIBLBFGS_SOURCES} ${BLOSSOM_SOURCES} ${WV_SOURCES} ${EGLITE_SOURCES} ${TET_SOURCES} ${MESHB_SOURCES} ${TRIANGLE_SOURCES} ${NUMPACK_SOURCES} ${JSON_SOURCES} ${GLAD_SOURCES} ${IMGUI_SOURCES} ${NNSEARCH_SOURCES} ${HLBFGS_SOURCES} ${QHULL_SOURCES} )
