ADD_SUBDIRECTORY( BlockLinAlg )
ADD_SUBDIRECTORY( DenseLinAlg )
ADD_SUBDIRECTORY( SparseLinAlg )

ADD_CUSTOM_TARGET( LinearAlgebraLib )
ADD_DEPENDENCIES( BlockLinAlgLib  LinearAlgebraLib )
ADD_DEPENDENCIES( DenseLinAlgLib  LinearAlgebraLib )
ADD_DEPENDENCIES( SparseLinAlgLib LinearAlgebraLib )
