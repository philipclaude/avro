#ifndef URSA_COMMON_TOOLS_H_
#define URSA_COMMON_TOOLS_H_

#include "common/types.h"

namespace ursa
{

#define UNUSED(x) (void)(x);

real_t randomWithin( const real_t lo , const real_t hi );
int randomWithin( const int lo , const int hi );

}

#endif
