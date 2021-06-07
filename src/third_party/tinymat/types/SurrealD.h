// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SURREALD_H
#define SURREALD_H

//#include <boost/type_traits/is_arithmetic.hpp>
#include <type_traits>

#include "PromoteSurreal.h"

// philip
#define SURREAL_TRAD

//#ifdef SURREAL_RVO
//#include "SurrealD_RVO.h"
#if defined(SURREAL_TRAD)
#include "SurrealD_Trad.h"
#elif defined(SURREAL_LAZY) || defined(SURREAL_REVERSE)
#include "SurrealD_Lazy.h"
#else
#error "Please define SURREAL_TRAD, SURREAL_LAZY or SURREAL_REVERSE"
#endif

//Let Surreal be treated as a arithmetic type
/*
philip
namespace boost
{
template<>
struct is_arithmetic<SurrealD> : boost::true_type {};
}
*/

namespace std
{
template<>
struct is_arithmetic<SurrealD> : std::true_type {};
}

#endif // SURREALD_H
