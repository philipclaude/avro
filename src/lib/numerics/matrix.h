//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_MATRIX_H_
#define avro_LIB_NUMERICS_MATRIX_H_

#include <vector>

//#include <tinymat/dense/dynamic/MatrixD.h>
//#include <tinymat/dense/dynamic/MatrixSymD.h>
//#include <tinymat/dense/dynamic/VectorD.h>

//#include <tinymat/dense/static/MatrixS.h>
//#include <tinymat/dense/static/VectorS.h>

//#include <tinymat/dense/tools/Identity.h>

#include "numerics/mat.h"
#include "numerics/symatd.h"

namespace avro
{

namespace numerics
{

#if 0
template<typename T> using MatrixD = tinymat::DLA::MatrixD<T>;
template<typename T> using SymMatrixD = tinymat::DLA::MatrixSymD<T>;

template<int M, int N, typename T> using MatrixS = tinymat::DLA::MatrixS<M,N,T>;
template<int M, typename T> using SymMatrixS = tinymat::DLA::MatrixSymS<M,T>;

template<int M,typename T> using VectorS = tinymat::DLA::VectorS<M,T>;
template<typename T> using VectorD = tinymat::DLA::VectorD<T>;
#else

template<typename T> using MatrixD = matd<T>;
template<typename T> using SymMatrixD = symd<T>;

template<int M, int N, typename T> using MatrixS = mats<M,N,T>;
//template<int M, typename T> using SymMatrixS = tinymat::DLA::MatrixSymS<M,T>;

template<int M,typename T> using VectorS = vecs<M,T>;
template<typename T> using VectorD = vecd<T>;

#endif

} // numerics

} // avro

#endif
