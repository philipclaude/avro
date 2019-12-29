if (luma_WITH_FPG)
	add_subdirectory(PCK/fpg)
endif()

add_source_directories( TET_SOURCES "tetgen source" tetgen1.5.0 )
set_source_files_properties( ${TET_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
set_source_files_properties( ${TET_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" TETLIBRARY )

add_definitions( -DTRILIBRARY )
add_source_directories( TRIANGLE_SOURCES "triangle source" triangle )
set_source_files_properties( ${TRIANGLE_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )

add_source_directories( TE_SOURCES "tiny expr source" tinyexpr )
set_source_files_properties( ${TE_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
include_directories( tinyexpr )

add_source_directories( JSON_SOURCES "json source" json )
set_source_files_properties( ${JSON_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
include_directories( json )

if (luma_BUILTIN_GLFW)
	add_subdirectory( glfw )
	set( THIRD_PARTY_INCLUDE_DIRS ${THIRD_PARTY_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/third_party/glfw/include )
	include_directories( ${CMAKE_SOURCE_DIR}/src/third_party/glfw/include )
endif()

add_source_directories( IMGUI_SOURCES "imgui source" imgui )
add_source_directories( IMGUI_SOURCES "imgui source" imgui/GL )
set_source_files_properties( ${IMGUI_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" IMGUI_IMPL_OPENGL_LOADER_GLAD )
include_directories( imgui )

set( GLAD_SOURCES )
if (luma_WITH_GL)
	add_source_directories( GLAD_SOURCES "glad source" glad/src )
	set( THIRD_PARTY_INCLUDE_DIRS ${THIRD_PARTY_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/src/third_party/glad/include )
	include_directories( ${CMAKE_SOURCE_DIR}/src/third_party/glad/include )
endif()

add_source_directories( MESHB_SOURCES "gmf source" libmeshb )
set_source_files_properties( ${MESHB_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
include_directories( libmeshb )

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

add_library(luma_third_party OBJECT ${SOURCES} ${TET_SOURCES} ${MESHB_SOURCES} ${TE_SOURCES} ${TRIANGLE_SOURCES} ${NUMPACK_SOURCES} ${JSON_SOURCES} ${GLAD_SOURCES} ${IMGUI_SOURCES} )
