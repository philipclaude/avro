#ifndef URSA_LIB_NUMERICS_MATRIX_H_
#define URSA_LIB_NUMERICS_MATRIX_H_

#include <vector>

#include <numpack/DenseLinAlg/DynamicSize/MatrixD.h>

namespace ursa
{

namespace numerics
{

template<typename T> using Matrix = numpack::DLA::MatrixD<T>;

} // numerics

} // ursa

#endif
