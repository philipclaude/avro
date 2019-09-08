// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MINMAX_H
#define MINMAX_H

//Some compilers make min/max macros...
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#include <type_traits>
//#include <boost/type_traits/is_arithmetic.hpp>

template<typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, T>::type
min(const T a, const T b) { return MIN(a,b); }

template<typename T>
inline typename std::enable_if< std::is_arithmetic<T>::value, T>::type
max(const T a, const T b) { return MAX(a,b); }

#endif //MINMAX_H
