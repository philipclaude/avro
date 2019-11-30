#ifndef LUNA_LIB_LINEAR_ALGEBRA_H_
#define LUNA_LIB_LINEAR_ALGEBRA_H_

#include "numerics/matrix.h"

#include <numpack/dense/dynamic/MatrixD_Det.h>

namespace luna
{

namespace numerics
{

template<typename type>
int
kernel( const MatrixD<type>& A , MatrixD<type>& K );

template<typename type>
type
determinant( const MatrixD<type>& A )
{
  return numpack::DLA::Det(A);
}

} // numerics

} // luna

#endif
