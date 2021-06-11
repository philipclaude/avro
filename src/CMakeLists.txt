# reset the warnings that were saved in third_party/CMakeLists.txt
set(CMAKE_C_FLAGS ${C_FLAGS0})
set(CMAKE_CXX_FLAGS ${CXX_FLAGS0})

# add all the source files
add_source_directories(SOURCES "Source Files" .)
add_source_directories(SOURCES "Source Files\\api" api)
add_source_directories(SOURCES "Source Files\\adaptation" lib/adaptation)
add_source_directories(SOURCES "Source Files\\common" lib/common)
add_source_directories(SOURCES "Source Files\\mesh" lib/element)
add_source_directories(SOURCES "Source Files\\geometry" lib/geometry)
add_source_directories(SOURCES "Source Files\\geometry-egads" lib/geometry/egads)
add_source_directories(SOURCES "Source Files\\geometry-psc" lib/geometry/psc)
add_source_directories(SOURCES "Source Files\\graphics" lib/graphics)
add_source_directories(SOURCES "Source Files\\library" lib/library)
add_source_directories(SOURCES "Source Files\\math" lib/math)
add_source_directories(SOURCES "Source Files\\mesh" lib/mesh)
add_source_directories(SOURCES "Source Files\\numerics" lib/numerics )
add_source_directories(SOURCES "Source Files\\numerics\\predicates" lib/numerics/predicates)
add_source_directories(SOURCES "Source Files\\voronoi" lib/voronoi)

# add the bin sources but remove the avro.cpp (main) executable source
add_source_directories(SOURCES "Source Files\\bin" bin)
list(REMOVE_ITEM SOURCES bin/avro.cpp )

#add_source_directories(SOURCES "Source Files\\api/python" api/python)

# create an object library for all the sources
add_library( avro_src OBJECT ${SOURCES} )

# create both shared and static libraries with same name
add_library(avro_shared SHARED $<TARGET_OBJECTS:avro_src> $<TARGET_OBJECTS:avro_third_party> )
set_target_properties( avro_shared PROPERTIES OUTPUT_NAME avro CLEAN_DIRECT_OUTPUT 1 )
add_library(avro_static STATIC $<TARGET_OBJECTS:avro_src> $<TARGET_OBJECTS:avro_third_party>)
set_target_properties( avro_static PROPERTIES OUTPUT_NAME avro CLEAN_DIRECT_OUTPUT 1 )

# create a target to build both shared and static libraries
add_custom_target( avro_lib DEPENDS avro_shared avro_static )

# add include directories for third party softwares
include_directories( ${ESP_INCLUDE_DIRS} )
include_directories( ${NLOPT_INCLUDE_DIRS} )
include_directories( ${THIRD_PARTY_INCLUDE_DIRS} )

include_directories( third_party/glad/include )
include_directories( third_party/glfw/include )
include_directories( third_party/tinymat )

# add the symbolic links for the webplot application
#install_symlink( ${ESP_WV_DIR} ${CMAKE_SOURCE_DIR}/app/WebViewer )

if (AVRO_NO_ESP)
  set( ESP_LIBRARIES )
  set (OCC_LIBRARIES )
  include_directories( third_party/egads )
  include_directories( third_party/wv )
  add_definitions( -DAVRO_NO_ESP )
endif()

set(avro_EXTERNAL_LIBRARIES ${ESP_LIBRARIES} ${OCC_LIBRARIES} ${LAPACK_LIBRARIES} ${NLOPT_LIBRARIES} ${MPI_LIBRARIES} ${PARMETIS_LIBRARIES} ${SUITESPARSE_LIBRARIES} dl ${BLAS_LIBRARIES} z )
if (avro_WITH_OPENMP)
  set( avro_EXTERNAL_LIBRARIES ${avro_EXTERNAL_LIBRARIES} ${OPENMP_LINK_FLAG} )
endif()

if (avro_MPI_ON)
  set( avro_EXTERNAL_LIBRARIES ${avro_EXTERNAL_LIBRARIES}
                               ${PARMETIS_LIBRARIES}
                               ${METIS_LIBRARIES}
                               ${MPI_CXX_LIBRARIES} )
  include_directories( ${PARMETIS_INCLUDE_DIRS} )
endif()

if (avro_OPENCL_ON)
  set( avro_EXTERNAL_LIBRARIES ${avro_EXTERNAL_LIBRARIES} ${OPENCL_LIBRARY} )
endif()

if (avro_WITH_GL)
  set( avro_EXTERNAL_LIBRARIES ${avro_EXTERNAL_LIBRARIES} ${OPENGL_LIBRARIES} ${OPENGL_gl_LIBRARY} glfw )
  if (APPLE)
    set( avro_EXTERNAL_LIBRARIES ${avro_EXTERNAL_LIBRARIES} ${COCOA_LIBRARY} )
  endif()
endif()

set( avro_EXTERNAL_LIBRARIES_SHARED ${avro_EXTERNAL_LIBRARIES}  dl )
set( avro_EXTERNAL_LIBRARIES_STATIC ${avro_EXTERNAL_LIBRARIES} ${ESP_LIBRARIES_STATIC} )

# add links to third party softwares
if (avro_STATIC_LINK_CXX)
  target_link_libraries( avro_shared ${avro_EXTERNAL_LIBRARIES_SHARED} -static-libgcc -static-libstdc++ )
else()
  target_link_libraries( avro_shared ${avro_EXTERNAL_LIBRARIES_SHARED})
endif()
target_link_libraries( avro_static ${avro_EXTERNAL_LIBRARIES_STATIC} )
