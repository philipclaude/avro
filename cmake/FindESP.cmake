

find_library( ESP_LIBRARY egads HINTS $ENV{ESP_DIR}/lib REQUIRED )
if (NOT ESP_LIBRARY)
	message( FATAL_ERROR "could not find ESP, try setting ESP_DIR" )
endif()

find_library( OCC_LIBRARY TKernel HINTS $ENV{CAS_DIR}/lib REQUIRED)
if (NOT OCC_LIBRARY)
	message( FATAL_ERROR "could not find OpenCASCADE library TKernel, try setting CAS_DIR" )
endif()

find_library( OCSM_LIBRARY ocsm HINTS $ENV{ESP_DIR}/lib REQUIRED)
if (NOT OCSM_LIBRARY)
	message( FATAL_ERROR "could not find OpenCSM library ocsm, try setting ESP_DIR" )
endif()


find_path( ESP_INCLUDE_DIR wsserver.h egads.h HINTS $ENV{ESP_DIR}/include REQUIRED )
find_path( ESP_WV_DIR wv-render.js HINTS $ENV{ESP_DIR}/wvClient/WebViewer REQUIRED )

get_filename_component( ESP_LIBRARY_PATH ${ESP_LIBRARY} PATH )
get_filename_component( OCC_LIBRARY_PATH ${OCC_LIBRARY} PATH )

set( OCC_LIBRARIES ${OCC_LIBRARY} )
unset( OCC_LIBRARY CACHE )

set( OCC_LIBRARY_PATH ${OCC_LIBRARY_PATH} CACHE PATH "OpenCASCADE library directory" )

set( ESP_LIBRARIES ${OCSM_LIBRARY} ${ESP_LIBRARY} wsserver emp  )
set( ESP_INCLUDE_DIRS ${ESP_INCLUDE_DIR} )

set( ESP_LIBRARIES_STATIC ${ESP_LIBRARY_PATH}/libegadstatic.a ${ESP_LIBRARY_PATH}/libwsserver.a )

set_property( DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY LINK_DIRECTORIES ${ESP_LIBRARY_PATH} ${OCC_LIBRARY_PATH} )

include_directories( ${ESP_INCLUDE_DIRS} )
