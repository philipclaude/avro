ADD_SUBDIRECTORY( block )
ADD_SUBDIRECTORY( dense )
ADD_SUBDIRECTORY( sparse )

ADD_CUSTOM_TARGET( LinearAlgebraLib )
ADD_DEPENDENCIES( blockLib  LinearAlgebraLib )
ADD_DEPENDENCIES( denseLib  LinearAlgebraLib )
ADD_DEPENDENCIES( sparseLib LinearAlgebraLib )
