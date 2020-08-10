include_directories( ${CMAKE_SOURCE_DIR}/src/third_party/tinymat )
add_executable( avro avro.cpp )
target_link_libraries( avro avro_shared ${AVRO_EXTERNAL_LIBRARIES} )
