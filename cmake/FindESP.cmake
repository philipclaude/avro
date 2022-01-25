set( HAVE_ESP TRUE )

find_library( ESP_LIBRARY egads HINTS $ENV{ESP_DIR}/lib )
if (NOT ESP_LIBRARY)
	#message( FATAL_ERROR "could not find ESP, try setting ESP_DIR" )
	set( HAVE_ESP FALSE )
endif()

find_library( OCC_LIBRARY TKernel HINTS $ENV{CAS_DIR}/lib )
if (NOT OCC_LIBRARY)
	#message( FATAL_ERROR "could not find OpenCASCADE library TKernel, try setting CAS_DIR" )
	set( HAVE_ESP FALSE )
endif()

find_library( OCSM_LIBRARY ocsm HINTS $ENV{ESP_DIR}/lib )
if (NOT OCSM_LIBRARY)
	#message( FATAL_ERROR "could not find OpenCSM library ocsm, try setting ESP_DIR" )
	set( HAVE_ESP FALSE )
endif()

if (NOT HAVE_ESP)
	message( STATUS "could not find the EngineeringSketchPad" )
	set(AVRO_NO_ESP TRUE)
else()

	find_path( ESP_INCLUDE_DIR egads.h HINTS $ENV{ESP_DIR}/include REQUIRED )

	get_filename_component( ESP_LIBRARY_PATH ${ESP_LIBRARY} PATH )
	get_filename_component( OCC_LIBRARY_PATH ${OCC_LIBRARY} PATH )

	set( OCC_LIBRARIES ${OCC_LIBRARY} )
	unset( OCC_LIBRARY CACHE )

	set( OCC_LIBRARY_PATH ${OCC_LIBRARY_PATH} CACHE PATH "OpenCASCADE library directory" )

	set( ESP_LIBRARIES ${OCSM_LIBRARY} ${ESP_LIBRARY} emp  )
	set( ESP_INCLUDE_DIRS ${ESP_INCLUDE_DIR} )

	set( ESP_LIBRARIES_STATIC ${ESP_LIBRARY_PATH}/libegadstatic.a )

	set_property( DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY LINK_DIRECTORIES ${ESP_LIBRARY_PATH} ${OCC_LIBRARY_PATH} )

	include_directories( ${ESP_INCLUDE_DIRS} )
endif()
