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

#include <numpack/dense/dynamic/MatrixD.h>
#include <numpack/dense/dynamic/MatrixSymD.h>
#include <numpack/dense/dynamic/VectorD.h>

#include <numpack/dense/static/MatrixS.h>
#include <numpack/dense/static/VectorS.h>

#include <numpack/dense/tools/Identity.h>

namespace avro
{

namespace numerics
{

template<typename T> using MatrixD = numpack::DLA::MatrixD<T>;
template<typename T> using SymMatrixD = numpack::DLA::MatrixSymD<T>;

template<int M, int N, typename T> using MatrixS = numpack::DLA::MatrixS<M,N,T>;
template<int M, typename T> using SymMatrixS = numpack::DLA::MatrixSymS<M,T>;

template<int M,typename T> using VectorS = numpack::DLA::VectorS<M,T>;
template<typename T> using VectorD = numpack::DLA::VectorD<T>;

using Mat4f = numpack::DLA::MatrixS<4,4,float>;
using Vec4f = numpack::DLA::VectorS<4,float>;

} // numerics

} // avro

#endif
