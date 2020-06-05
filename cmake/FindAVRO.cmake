
find_library( avro_LIBRARY avro HINTS $ENV{avro_DIR}/build/release/lib REQUIRED )
if (NOT avro_LIBRARY)
  message( STATUS "could not find avro" )
endif()

find_path( avro_INCLUDE_DIR avro.h HINTS $ENV{avro_DIR}/api REQUIRED )
get_filename_component( avro_LIBRARY_PATH ${avro_LIBRARY} PATH )

set( avro_INCLUDE_DIRS ${avro_INCLUDE_DIR} $ENV{avro_DIR}/src )

message( STATUS "found avro library at ${avro_LIBRARY_PATH}" )
message( STATUS "found avro headers at ${avro_INCLUDE_DIRS}" )
