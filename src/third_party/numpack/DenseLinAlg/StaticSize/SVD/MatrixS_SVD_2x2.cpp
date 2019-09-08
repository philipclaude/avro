// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MatrixS_SVD.h"
#include "../VectorS.h"
#include "../MatrixSymS.h"

#include "tools/SANSnumerics.h"

#include "numpack/types/SurrealS.h"

#include <cmath>
#include <limits>

namespace numpack 
{
namespace DLA
{


/**
  * @class JacobiRotation
  * @brief Rotation given by a cosine-sine pair.
  *
  * This class represents a Jacobi or Givens rotation.
  * This is a 2D rotation in the plane \c J of angle \f$ \theta \f$ defined by
  * its cosine \c c and sine \c s as follows:
  * \f$ J = \left ( \begin{array}{cc} c & \overline s \\ -s  & \overline c \end{array} \right ) \f$
  *
  * You can apply the respective counter-clockwise rotation to a column vector
  * \c v by* applying its adjoint on the left: \f$ v = J^* v \f$
  * \endcode
  */

template<typename T>
class Jacobi_rotation_2x2
{
public:
  // Default constructor without any initialization.
  Jacobi_rotation_2x2() {}

  /// Construct a planar rotation from a cosine-sine pair (\a c, \c s).
  Jacobi_rotation_2x2(const T& c, const T& s) : _cos(c), _sin(s) {}

        T& cos()       { return _cos; }
        T& sin()       { return _sin; }
  const T& cos() const { return _cos; }
  const T& sin() const { return _sin; }

  // Concatenates two planar rotation
  Jacobi_rotation_2x2 operator*(const Jacobi_rotation_2x2& jr)
  {
    return Jacobi_rotation_2x2(_cos * jr._cos - _sin * jr._sin,
                               _cos * jr._sin + _sin * jr._cos);
  }

  // Returns the transposed transformation
  Jacobi_rotation_2x2 transpose() const { return Jacobi_rotation_2x2(_cos, -_sin); }

  // Makes \c *this as a Jacobi rotation \c J such that applying \a J on both
  // the right and left sides of the 2x2 selfadjoint matrix
  // \f$ B = \left ( \begin{array}{cc} \text{this}_{pp} & \text{this}_{pq} \\ (\text{this}_{pq})^* & \text{this}_{qq} \end{array} \right )\f$
  // yields a diagonal matrix \f$ A = J^* B J \f$
  inline bool make_jacobi(MatrixS<2,2,T>& m, int p, int q)
  {
    return make_jacobi( m(p, p), m(p, q), m(q, q) );
  }

  // Makes \c *this as a Jacobi rotation \a J such that applying \a J on both
  // the right and left sides of the selfadjoint 2x2 matrix
  // \f$ B = \left ( \begin{array}{cc} x & y \\ \overline y & z \end{array} \right )\f$
  // yields a diagonal matrix \f$ A = J^* B J \f$
  bool make_jacobi(const T& x, const T& y, const T& z);

  // Applies the rotation in the plane \a j to the rows \a p and \a q
  // of \c m, i.e., it computes B = J * B,
  // with \f$ B = \left ( \begin{array}{cc} \text{m.row}(p) \\ \text{m.row}(q) \end{array} \right ) \f$.
  inline void apply_on_the_left(MatrixS<2,2,T>& m, int p, int q)
  {
    for (int i = 0; i < 2; ++i)
    {
      T xi = m(p, i);
      T yi = m(q, i);
      m(p, i) =  cos() * xi + sin() * yi;
      m(q, i) = -sin() * xi + cos() * yi;
    }
  }

  // Applies the rotation in the plane \a j to the columns \a p and \a q of
  // \c m, i.e., it computes B = B * J with
  // \f$ B = \left ( \begin{array}{cc} \text{m.col}(p) & \text{m.col}(q) \end{array} \right ) \f$.
  inline void apply_on_the_right(MatrixS<2,2,T>& m, int p, int q)
  {
    for (int i = 0; i < 2; ++i)
    {
      T xi = m(i, p);
      T yi = m(i, q);
      m(i, p) =  cos() * xi - sin() * yi;
      m(i, q) =  sin() * xi + cos() * yi;
    }
  }

private:
  T _cos, _sin;
};

// -----------------------------------------------------------------------------

template<typename T>
bool Jacobi_rotation_2x2<T>::make_jacobi(const T& x, const T& y, const T& z)
{
  if (y == 0.0)
  {
    _cos = 1.0;
    _sin = 0.0;
    return false;
  }
  else
  {
    T tau = (x-z) / (2.0 * fabs(y));
    T w = sqrt(tau*tau + 1.0);
    T t;
    if (tau > 0.0)
    {
      t = 1.0/(tau + w);
    }
    else
    {
      t = 1.0/(tau - w);
    }
    T sign_t = t > 0.0 ? 1.0 : -1.0;
    T n = 1.0/sqrt( t*t + 1.0 );
    _sin = -sign_t * (y/fabs(y)) * fabs(t) * n;
    _cos = n;
    return true;
  }
}

// -----------------------------------------------------------------------------

template<typename T>
void real_2x2_jacobi_svd(
        const MatrixS<2,2,T>& matrix,
        Jacobi_rotation_2x2<T>* j_left,
        Jacobi_rotation_2x2<T>* j_right)
{
  MatrixS<2,2,T> m = 0;

  m(0, 0) = matrix(1, 1); m(0, 1) = matrix(1, 0);
  m(1, 0) = matrix(0, 1); m(1, 1) = matrix(0, 0);

  Jacobi_rotation_2x2<T> rot;

  T t = m(0, 0) + m(1, 1);
  T d = m(1, 0) - m(0, 1);

  if (t == 0.0)
  {
    rot.cos() = 0.0;
    rot.sin() = d > 0.0 ? 1.0 : -1.0;
  }
  else
  {
    T u = d/t;
    rot.cos() = 1.0/sqrt(1.0 + u*u);
    rot.sin() = rot.cos() * u;
  }

  rot.apply_on_the_left(m, 0, 1);

  j_right->make_jacobi(m, 0, 1);

  *j_left = rot * j_right->transpose();
}


template< int M, int N, class T >
void
SVD( const MatrixS<M,N,T>& A, MatrixS<M,M,T>& U, VectorS<MIN(M,N),T>& S, MatrixS<N,N,T>& VT )
{
  // currently we stop when we reach precision 2*epsilon as the last bit of
  // precision can require an unreasonable number of iterations,
  // only worsening the precision of U and V as we accumulate more rotations
  const Real precision = 2.0 * std::numeric_limits<Real>::epsilon();

  // limit for very small denormal numbers to be considered zero in order to
  // avoid infinite loops
  const Real consider_null = 2.0 * std::numeric_limits<Real>::denorm_min();

  MatrixS<M,N,T> work_mat = A;

  U = Identity();
  MatrixS<N,N,T> V = Identity();
  S = 0.0;

  // The main Jacobi SVD iteration
  // if this 2x2 matrix is not diagonal already...
  // notice that this comparison will evaluate to false if any NaN is involved,
  // similarly, small denormal numbers are considered zero.
  T max_diag = max( fabs(work_mat(1, 1)), fabs(work_mat(0, 0)) );
  T tmp = precision*max_diag;

  T threshold = max(consider_null, tmp);

  T max_anti_diag = max( fabs(work_mat(1, 0)), fabs(work_mat(0, 1)) );

  if ( max_anti_diag > threshold)
  {
    Jacobi_rotation_2x2<T> j_left, j_right;
    real_2x2_jacobi_svd(work_mat, &j_left, &j_right);

    // accumulate resulting Jacobi rotations
    j_left.apply_on_the_left(work_mat, 1, 0);
    (j_left.transpose()).apply_on_the_right(U, 1, 0);

    j_right.apply_on_the_right(work_mat, 1, 0);
    j_right.apply_on_the_right(V, 1, 0);
  }

  // The work matrix is now diagonal,
  // so ensure it's positive so its diagonal entries are the singular values
  for (int i = 0; i < 2; ++i)
  {
    T a = fabs( work_mat(i, i) );
    S[i] = a;

    if ( a != 0.0 )
    {
      VectorS<2,T> col = U.col(i) * work_mat(i, i)/a;
      U(0,i) = col[0];
      U(1,i) = col[1];
    }
  }

  // Sort singular values in descending order
  if ( S[0] < S[1] )
  {
    //Swap singular values
    T tmp = S[0];
    S[0] = S[1];
    S[1] = tmp;

    //Swap columns of unitary matrices
    U.swap_cols(0,1);
    V.swap_cols(0,1);
  }

  //Need to return the transpose of V
  VT = Transpose(V);
}


#define INSTANTIATE_SVD(T) \
template void SVD<2,2,T>(const MatrixS<2,2,T>& A, MatrixS<2,2,T>& U, VectorS<2,T>& S, MatrixS<2,2,T>& VT);

INSTANTIATE_SVD(Real)
INSTANTIATE_SVD(SurrealS<1>)

}
}
