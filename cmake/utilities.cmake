
# utility for adding an executable
macro( luna_add_executable )

    if( NOT LUNA_BUILD_DYNAMIC)
        #set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static" )
        #set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -static" )
    endif()
    add_executable( ${ARGN} )
    add_dependencies( ${ARGV0} luna_lib )
    target_link_libraries( ${ARGV0} luna ${LUNA_EXTERNAL_LIBRARIES} )
endmacro()

#!
# @brief Add sources from directories
# @details
# Add the sources from the specified \p directories to variable \p var
# and place them in the specified \p folder if non empty.
# @param[out] var name of the variable that receives the result list
# @param[in] folder the name of the folder
# @param[in] directories list of directories to scan
#
function(add_source_directories var folder)
    set(sources)
    foreach(dir ${ARGN})
        file(GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${dir}/*.[ch]" "${dir}/*.[ch]pp" "${dir}/*.[ch]xx")
        list(APPEND sources ${_sources})
    endforeach()

    if( NOT folder STREQUAL "")
        source_group(${folder} FILES ${sources})
    endif()

    set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

function(add_test_files var type )
  set(sources)
  foreach(dir ${ARGN})
      file(GLOB _sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "${dir}/*_${type}.cpp")
      list(APPEND sources ${_sources})
  endforeach()
  set(${var} ${${var}} ${sources} PARENT_SCOPE)
endfunction()

# utility for adding a symbolic link
macro(install_symlink filepath sympath)
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${filepath} ${sympath})
endmacro(install_symlink)

function( add_test testname testfile ${ARGN} )
  set( TEST_NAMES ${TEST_NAMES} ${testname} PARENT_SCOPE )
  set( TEST_FILES0 ${TEST_FILES} )
  list( LENGTH TEST_FILES0 FILESTART )
  set( TEST_FILES0 ${TEST_FILES} ${testfile} )
  list( LENGTH TEST_FILES0 FILEEND )
  set( TEST_FILES ${TEST_FILES0} PARENT_SCOPE )
  set(TEST_FILES_START ${TEST_FILES_START} ${FILESTART} PARENT_SCOPE )
  set(TEST_FILES_END ${TEST_FILES_END} ${FILEEND} PARENT_SCOPE )
  set(numprocs ${ARGN})
  list( LENGTH numprocs num_extra_args )
  if ( ${num_extra_args} GREATER 0 )
    set(TEST_NUMPROCS ${TEST_NUMPROCS} ${numprocs} PARENT_SCOPE )
  else()
    set(TEST_NUMPROCS ${TEST_NUMPROCS} 0 PARENT_SCOPE )
  endif()
endfunction()

function( get_testname i testname )
  # loop through each testname
  foreach( name ${TEST_NAMES} )

    # get the bounding indices for this testname
    list( FIND TEST_NAMES ${name} j )
    list( GET TEST_FILES_START ${j} f0 )
    list( GET TEST_FILES_END ${j} f1 )

    if ( ${f0} EQUAL ${i} )
      set( ${testname} ${name} PARENT_SCOPE )
      return()
    endif()

    if ( ${f0} LESS ${i} )
      if (${f1} GREATER ${i} )
        set( ${testname} ${name} PARENT_SCOPE )
        return()
      endif()
    endif()
  endforeach()
endfunction()

function( get_testname_fromfile filename testname )
  list( FIND TEST_FILES ${filename} i )
  get_testname( ${i} testname0 )
  set( ${testname} ${testname0} PARENT_SCOPE )
endfunction()
