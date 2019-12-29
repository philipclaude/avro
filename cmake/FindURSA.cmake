
find_library( luma_LIBRARY luma HINTS $ENV{luma_DIR}/build/release/lib REQUIRED )
if (NOT luma_LIBRARY)
  message( STATUS "could not find luma" )
endif()

find_path( luma_INCLUDE_DIR luma.h HINTS $ENV{luma_DIR}/api REQUIRED )
get_filename_component( luma_LIBRARY_PATH ${luma_LIBRARY} PATH )

set( luma_INCLUDE_DIRS ${luma_INCLUDE_DIR} $ENV{luma_DIR}/src )

message( STATUS "found luma library at ${luma_LIBRARY_PATH}" )
message( STATUS "found luma headers at ${luma_INCLUDE_DIRS}" )
