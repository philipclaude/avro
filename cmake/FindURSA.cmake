
find_library( URSA_LIBRARY ursa HINTS $ENV{URSA_DIR}/build/release/lib REQUIRED )
if (NOT URSA_LIBRARY)
  message( STATUS "could not find ursa" )
endif()

find_path( URSA_INCLUDE_DIR ursa.h HINTS $ENV{URSA_DIR}/api REQUIRED )
get_filename_component( URSA_LIBRARY_PATH ${URSA_LIBRARY} PATH )

set( URSA_INCLUDE_DIRS ${URSA_INCLUDE_DIR} $ENV{URSA_DIR}/src )

message( STATUS "found ursa library at ${URSA_LIBRARY_PATH}" )
message( STATUS "found ursa headers at ${URSA_INCLUDE_DIRS}" )
