// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_NUMERICS_DETERMINANT_H_
#define AVRO_NUMERICS_DETERMINANT_H_

#include "common/types.h"

namespace avro
{

template<typename type> class matd;

namespace numerics
{

template <typename T> class symd;

template<typename type> type det(const matd<type>& X);
template<typename type> type det(const symd<type>& X);

} // numerics

} // avro

#endif
