//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_DETERMINANT_H_
#define avro_LIB_NUMERICS_DETERMINANT_H_

#include "common/types.h"

namespace avro
{

namespace numerics
{

template<typename type> class densMat;

template<typename type> type determinant(const densMat<type>& X);

} // numerics

} // avro

#endif
