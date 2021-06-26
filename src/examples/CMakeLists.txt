add_executable( example_adapt_serial adapt_serial.cpp )
target_link_libraries( example_adapt_serial avro_shared ${AVRO_EXTERNAL_LIBRARIES} )

add_executable( example_adapt_parallel adapt_parallel.cpp )
target_link_libraries( example_adapt_parallel avro_shared ${AVRO_EXTERNAL_LIBRARIES} )
