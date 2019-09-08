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

#add_source_directories( GMF_SOURCES "gmf source" libMeshb/libMeshb7 )
#set_source_files_properties( ${GMF_SOURCES} PROPERTIES "COMPILE_DEFINITIONS" TRANSMESH )
#set_source_files_properties( ${GMF_SOURCES} PROPERTIES COMPILE_FLAGS "-w" )
#include_directories( libMeshb/libMeshb7 )

add_library(ursa_third_party OBJECT ${SOURCES} ${TET_SOURCES} ${ANN_SOURCES} ${GMF_SOURCES} ${TE_SOURCES} ${TRIANGLE_SOURCES} )
