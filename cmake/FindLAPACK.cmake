
if( APPLE )
    # OSX provides the lapack library in the Accelerate framewwork
    find_library(LAPACK_LIBRARIES NAMES Accelerate )
    find_path(LAPACKE_INCLUDE_DIRS Accelerate/Accelerate.h)
    add_definitions( -DLAPACK_ACCELERATE )
else() 
    find_library( LAPACK_LIBRARIES NAMES lapack )
endif()
