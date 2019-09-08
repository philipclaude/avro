// ursa: unstructured adaptation library 
// Copyright 2019-present, Philip Claude Caplan
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef URSA_COMMON_TOOLS_H_
#define URSA_COMMON_TOOLS_H_

#include "common/types.h"

namespace ursa
{

#define UNUSED(x) (void)(x);

real randomWithin( const real lo , const real hi );
int randomWithin( const int lo , const int hi );

}

#endif
