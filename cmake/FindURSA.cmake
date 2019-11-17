
find_library( LUNA_LIBRARY luna HINTS $ENV{LUNA_DIR}/build/release/lib REQUIRED )
if (NOT LUNA_LIBRARY)
  message( STATUS "could not find luna" )
endif()

find_path( LUNA_INCLUDE_DIR luna.h HINTS $ENV{LUNA_DIR}/api REQUIRED )
get_filename_component( LUNA_LIBRARY_PATH ${LUNA_LIBRARY} PATH )

set( LUNA_INCLUDE_DIRS ${LUNA_INCLUDE_DIR} $ENV{LUNA_DIR}/src )

message( STATUS "found luna library at ${LUNA_LIBRARY_PATH}" )
message( STATUS "found luna headers at ${LUNA_INCLUDE_DIRS}" )
