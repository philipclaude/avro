add_executable( avro avro.cpp )
target_link_libraries( avro avro_shared ${AVRO_EXTERNAL_LIBRARIES} )
