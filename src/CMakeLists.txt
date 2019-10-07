# reset the warnings that were saved in third_party/CMakeLists.txt
set(CMAKE_C_FLAGS ${C_FLAGS0})
set(CMAKE_CXX_FLAGS ${CXX_FLAGS0})

# add all the source files
add_source_directories(SOURCES "Source Files" .)
add_source_directories(SOURCES "Source Files\\common" lib/common)
add_source_directories(SOURCES "Source Files\\geometrics" lib/geometrics)
add_source_directories(SOURCES "Source Files\\graphics" lib/graphics)
add_source_directories(SOURCES "Source Files\\library" lib/library)
add_source_directories(SOURCES "Source Files\\mesh" lib/master)
add_source_directories(SOURCES "Source Files\\mesh" lib/mesh)
add_source_directories(SOURCES "Source Files\\numerics" lib/numerics )
add_source_directories(SOURCES "Source Files\\numerics\\predicates" lib/numerics/predicates)


# add the bin sources but remove the ursa.cpp (main) executable source
add_source_directories(SOURCES "Source Files\\bin" bin)
list(REMOVE_ITEM SOURCES bin/ursa.cpp )

# create an object library for all the sources
add_library( ursa_src OBJECT ${SOURCES} )

# create both shared and static libraries with same name
add_library(ursa_shared SHARED $<TARGET_OBJECTS:ursa_src> $<TARGET_OBJECTS:ursa_third_party> )
set_target_properties( ursa_shared PROPERTIES OUTPUT_NAME ursa CLEAN_DIRECT_OUTPUT 1 )
add_library(ursa_static STATIC $<TARGET_OBJECTS:ursa_src> $<TARGET_OBJECTS:ursa_third_party>)
set_target_properties( ursa_static PROPERTIES OUTPUT_NAME ursa CLEAN_DIRECT_OUTPUT 1 )

# create a target to build both shared and static libraries
add_custom_target( ursa_lib DEPENDS ursa_shared ursa_static )

# add include directories for third party softwares
include_directories( ${ESP_INCLUDE_DIRS} )
include_directories( ${NLOPT_INCLUDE_DIRS} )
include_directories( ${THIRD_PARTY_INCLUDE_DIRS} )
message( STATUS "tp incldue = " ${THIRD_PARTY_INCLUDE_DIRS} )

include_directories( third_party/glad/include )
include_directories( third_party/glfw/include )

# add the symbolic links for the webplot application
#install_symlink( ${ESP_WV_DIR} ${CMAKE_SOURCE_DIR}/app/WebViewer )

set(URSA_EXTERNAL_LIBRARIES ${ESP_LIBRARIES} ${OCC_LIBRARIES} ${LAPACK_LIBRARIES} ${NLOPT_LIBRARIES} dl blas z )
if (URSA_WITH_OPENMP)
  set( URSA_EXTERNAL_LIBRARIES ${URSA_EXTERNAL_LIBRARIES} ${OPENMP_LINK_FLAG} )
endif()

if (URSA_MPI_ON)
  set( URSA_EXTERNAL_LIBRARIES ${URSA_EXTERNAL_LIBRARIES}
                               ${PARMETIS_LIBRARIES}
                               ${METIS_LIBRARIES}
                               ${MPI_CXX_LIBRARIES} )
  include_directories( ${PARMETIS_INCLUDE_DIRS} )
endif()

if (URSA_OPENCL_ON)
  set( URSA_EXTERNAL_LIBRARIES ${URSA_EXTERNAL_LIBRARIES} ${OPENCL_LIBRARY} )
endif()

if (URSA_WITH_GL)
  set( URSA_EXTERNAL_LIBRARIES ${URSA_EXTERNAL_LIBRARIES} ${OPENGL_LIBRARIES} ${OPENGL_gl_LIBRARY} glfw )
  if (APPLE)
    set( URSA_EXTERNAL_LIBRARIES ${URSA_EXTERNAL_LIBRARIES} ${COCOA_LIBRARY} )
  endif()
endif()

set( URSA_EXTERNAL_LIBRARIES_SHARED ${URSA_EXTERNAL_LIBRARIES}  dl )
set( URSA_EXTERNAL_LIBRARIES_STATIC ${URSA_EXTERNAL_LIBRARIES} ${ESP_LIBRARIES_STATIC} )

# add links to third party softwares
if (URSA_STATIC_LINK_CXX)
  target_link_libraries( ursa_shared ${URSA_EXTERNAL_LIBRARIES_SHARED} -static-libgcc -static-libstdc++ )
else()
  target_link_libraries( ursa_shared ${URSA_EXTERNAL_LIBRARIES_SHARED})
endif()
target_link_libraries( ursa_static ${URSA_EXTERNAL_LIBRARIES_STATIC} )
