#ifndef URSA_LIB_NUMERICS_MATRIX_H_
#define URSA_LIB_NUMERICS_MATRIX_H_

#include <vector>

namespace ursa
{

namespace numerics
{

template<typename T> class Vector;

template<typename T>
class MatrixBase
{

protected:
  std::vector<T> data_;
};

template<typename T>
class DenseMatrix : public MatrixBase<T>
{

};

template<typename T>
class SparseMatrix
{

};

template<typename T>
class SymmetricMatrix
{

};

} // numerics

} // ursa

#endif
