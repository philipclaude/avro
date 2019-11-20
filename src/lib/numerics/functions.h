#ifndef LUNA_NUMERICS_FUNCTIONS_H_
#define LUNA_NUMERICS_FUNCTIONS_H_

#include "common/types.h"

namespace luna
{

namespace numerics
{

index_t factorial( const index_t i );
index_t binomial( const index_t n , const index_t k );
index_t nchoosek( index_t n , index_t k );

} // numerics

} // luna

#endif
