#ifndef LUNA_LIB_NUMERICS_MATRIX_H_
#define LUNA_LIB_NUMERICS_MATRIX_H_

#include <vector>

#include <numpack/dense/dynamic/MatrixD.h>
#include <numpack/dense/dynamic/MatrixSymD.h>
#include <numpack/dense/dynamic/VectorD.h>

#include <numpack/dense/static/MatrixS.h>
#include <numpack/dense/static/VectorS.h>

#include <numpack/dense/tools/Identity.h>

namespace luna
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

} // luna

#endif
