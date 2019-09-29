#ifndef URSA_LIB_NUMERICS_MATRIX_H_
#define URSA_LIB_NUMERICS_MATRIX_H_

#include <vector>

#include <numpack/dense/dynamic/MatrixD.h>
#include <numpack/dense/dynamic/MatrixSymD.h>

namespace ursa
{

namespace numerics
{

template<typename T> using MatrixD = numpack::DLA::MatrixD<T>;
template<typename T> using SymMatrixD = numpack::DLA::MatrixSymD<T>;

template<int M, int N, typename T> using MatrixS = numpack::DLA::MatrixS<M,N,T>;
template<int M, typename T> using SymMatrixS = numpack::DLA::MatrixSymS<M,T>;

using Mat4f = numpack::DLA::MatrixS<4,4,float>;
using Vec4f = numpack::DLA::VectorS<4,float>;

} // numerics

} // ursa

#endif
