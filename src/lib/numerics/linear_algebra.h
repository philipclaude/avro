#ifndef luma_LIB_LINEAR_ALGEBRA_H_
#define luma_LIB_LINEAR_ALGEBRA_H_

#include "common/error.h"

#include "numerics/matrix.h"

#include <numpack/dense/dynamic/MatrixD_Det.h>
#include <numpack/dense/dynamic/MatrixD_Diag.h>
#include <numpack/dense/dynamic/Eigen.h>
#include <numpack/dense/InverseLU.h>
#include <numpack/Transpose.h>

#include <numpack/types/SurrealD.h>
#include <numpack/types/SurrealS.h>

#include <cmath>

namespace luma
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

template<typename type>
type
determinant( const SymMatrixD<type>& A )
{
  return numpack::DLA::Det(A);
}

template< class T >
inline numpack::DLA::MatrixSymD<T>
logm(const numpack::DLA::MatrixSymD<T>& A)
{
  // compute the eigensystem with normalized eigenvectors
  numpack::DLA::EigenSystemPair<T> LE(A);

  // compute the log of the eigenvalues
  for (int i = 0; i < A.m(); i++ )
    LE.L[i] = ::log(LE.L[i]);

  // return the symmetric matrix
  return LE;
}

template< class T >
inline numpack::DLA::MatrixSymD<T>
expm(const numpack::DLA::MatrixSymD<T>& A)
{
  // compute the eigensystem with normalized eigenvectors
  numpack::DLA::EigenSystemPair<T> LE(A);

  // compute the log of the eigenvalues
  for (int i = 0; i < A.m(); i++ )
    LE.L[i] = ::exp(LE.L[i]);

  // return the symmetric matrix
  return LE;
}

template< class T >
inline numpack::DLA::MatrixSymD<T>
sqrtm(const numpack::DLA::MatrixSymD<T>& A)
{
  // compute the eigensystem with normalized eigenvectors
  numpack::DLA::EigenSystemPair<T> LE(A);

  // compute the log of the eigenvalues
  for (int i = 0; i < A.m(); i++ )
    LE.L[i] = ::sqrt(LE.L[i]);

  // return the symmetric matrix
  return LE;
}

template<typename T>
inline MatrixD<T>
inverse( MatrixD<T>& M )
{
	luma_assert( M.m() == M.n() );
	T idetM = 1./determinant(M);

  MatrixD<T> Minv(M.m(),M.n());

	if (M.n()==1)
	{
		Minv(0,0) = idetM;
	}
	else if (M.n()==2)
	{
		Minv(0,0) =  M(1,1)*idetM;
		Minv(0,1) = -M(0,1)*idetM;
		Minv(1,0) = -M(1,0)*idetM;
		Minv(1,1) =  M(0,0)*idetM;
	}
	else if (M.n()==3)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2);
		Minv(0,0) = (a2_2*a3_3 -a2_3*a3_2)*idetM;
		Minv(0,1) = (a1_3*a3_2 -a1_2*a3_3)*idetM;
		Minv(0,2) = (a1_2*a2_3 -a1_3*a2_2)*idetM;
		Minv(1,0) = (a2_3*a3_1 -a2_1*a3_3)*idetM;
		Minv(1,1) = (a1_1*a3_3 -a1_3*a3_1)*idetM;
		Minv(1,2) = (a1_3*a2_1 -a1_1*a2_3)*idetM;
		Minv(2,0) = (a2_1*a3_2 -a2_2*a3_1)*idetM;
		Minv(2,1) = (a1_2*a3_1 -a1_1*a3_2)*idetM;
		Minv(2,2) = (a1_1*a2_2 -a1_2*a2_1)*idetM;
	}
	else if (M.n()==4)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2); T a1_4 = M(0,3);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2); T a2_4 = M(1,3);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2); T a3_4 = M(2,3);
		T a4_1 = M(3,0); T a4_2 = M(3,1); T a4_3 = M(3,2); T a4_4 = M(3,3);

		Minv(0,0) = (a2_2*a3_3*a4_4-a2_2*a3_4*a4_3-a2_3*a3_2*a4_4+a2_3*a3_4*a4_2+a2_4*a3_2*a4_3-a2_4*a3_3*a4_2)*idetM;
		Minv(0,1) = (-a1_2*a3_3*a4_4+a1_2*a3_4*a4_3+a1_3*a3_2*a4_4-a1_3*a3_4*a4_2-a1_4*a3_2*a4_3+a1_4*a3_3*a4_2)*idetM;
		Minv(0,2) = (a1_2*a2_3*a4_4-a1_2*a2_4*a4_3-a1_3*a2_2*a4_4+a1_3*a2_4*a4_2+a1_4*a2_2*a4_3-a1_4*a2_3*a4_2)*idetM;
		Minv(0,3) = (-a1_2*a2_3*a3_4+a1_2*a2_4*a3_3+a1_3*a2_2*a3_4-a1_3*a2_4*a3_2-a1_4*a2_2*a3_3+a1_4*a2_3*a3_2)*idetM;
		Minv(1,0) = (-a2_1*a3_3*a4_4+a2_1*a3_4*a4_3+a2_3*a3_1*a4_4-a2_3*a3_4*a4_1-a2_4*a3_1*a4_3+a2_4*a3_3*a4_1)*idetM;
		Minv(1,1) = (a1_1*a3_3*a4_4-a1_1*a3_4*a4_3-a1_3*a3_1*a4_4+a1_3*a3_4*a4_1+a1_4*a3_1*a4_3-a1_4*a3_3*a4_1)*idetM;
		Minv(1,2) = (-a1_1*a2_3*a4_4+a1_1*a2_4*a4_3+a1_3*a2_1*a4_4-a1_3*a2_4*a4_1-a1_4*a2_1*a4_3+a1_4*a2_3*a4_1)*idetM;
		Minv(1,3) = (a1_1*a2_3*a3_4-a1_1*a2_4*a3_3-a1_3*a2_1*a3_4+a1_3*a2_4*a3_1+a1_4*a2_1*a3_3-a1_4*a2_3*a3_1)*idetM;
		Minv(2,0) = (a2_1*a3_2*a4_4-a2_1*a3_4*a4_2-a2_2*a3_1*a4_4+a2_2*a3_4*a4_1+a2_4*a3_1*a4_2-a2_4*a3_2*a4_1)*idetM;
		Minv(2,1) = (-a1_1*a3_2*a4_4+a1_1*a3_4*a4_2+a1_2*a3_1*a4_4-a1_2*a3_4*a4_1-a1_4*a3_1*a4_2+a1_4*a3_2*a4_1)*idetM;
		Minv(2,2) = (a1_1*a2_2*a4_4-a1_1*a2_4*a4_2-a1_2*a2_1*a4_4+a1_2*a2_4*a4_1+a1_4*a2_1*a4_2-a1_4*a2_2*a4_1)*idetM;
		Minv(2,3) = (-a1_1*a2_2*a3_4+a1_1*a2_4*a3_2+a1_2*a2_1*a3_4-a1_2*a2_4*a3_1-a1_4*a2_1*a3_2+a1_4*a2_2*a3_1)*idetM;
		Minv(3,0) = (-a2_1*a3_2*a4_3+a2_1*a3_3*a4_2+a2_2*a3_1*a4_3-a2_2*a3_3*a4_1-a2_3*a3_1*a4_2+a2_3*a3_2*a4_1)*idetM;
		Minv(3,1) = (a1_1*a3_2*a4_3-a1_1*a3_3*a4_2-a1_2*a3_1*a4_3+a1_2*a3_3*a4_1+a1_3*a3_1*a4_2-a1_3*a3_2*a4_1)*idetM;
		Minv(3,2) = (-a1_1*a2_2*a4_3+a1_1*a2_3*a4_2+a1_2*a2_1*a4_3-a1_2*a2_3*a4_1-a1_3*a2_1*a4_2+a1_3*a2_2*a4_1)*idetM;
		Minv(3,3) = (a1_1*a2_2*a3_3-a1_1*a2_3*a3_2-a1_2*a2_1*a3_3+a1_2*a2_3*a3_1+a1_3*a2_1*a3_2-a1_3*a2_2*a3_1)*idetM;
	}
	else
		luma_implement;
  return Minv;
}

template<typename T>
inline SymMatrixD<T>
inverse( SymMatrixD<T>& M )
{
	T idetM = 1./determinant(M);

  SymMatrixD<T> Minv(M.m(),M.n());

	if (M.n()==1)
	{
		Minv(0,0) = idetM;
	}
	else if (M.n()==2)
	{
		Minv(0,0) =  M(1,1)*idetM;
		Minv(0,1) = -M(0,1)*idetM;
		Minv(1,0) = -M(1,0)*idetM;
		Minv(1,1) =  M(0,0)*idetM;
	}
	else if (M.n()==3)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2);
		Minv(0,0) = (a2_2*a3_3 -a2_3*a3_2)*idetM;
		Minv(0,1) = (a1_3*a3_2 -a1_2*a3_3)*idetM;
		Minv(0,2) = (a1_2*a2_3 -a1_3*a2_2)*idetM;
		Minv(1,0) = (a2_3*a3_1 -a2_1*a3_3)*idetM;
		Minv(1,1) = (a1_1*a3_3 -a1_3*a3_1)*idetM;
		Minv(1,2) = (a1_3*a2_1 -a1_1*a2_3)*idetM;
		Minv(2,0) = (a2_1*a3_2 -a2_2*a3_1)*idetM;
		Minv(2,1) = (a1_2*a3_1 -a1_1*a3_2)*idetM;
		Minv(2,2) = (a1_1*a2_2 -a1_2*a2_1)*idetM;
	}
	else if (M.n()==4)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2); T a1_4 = M(0,3);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2); T a2_4 = M(1,3);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2); T a3_4 = M(2,3);
		T a4_1 = M(3,0); T a4_2 = M(3,1); T a4_3 = M(3,2); T a4_4 = M(3,3);

		Minv(0,0) = (a2_2*a3_3*a4_4-a2_2*a3_4*a4_3-a2_3*a3_2*a4_4+a2_3*a3_4*a4_2+a2_4*a3_2*a4_3-a2_4*a3_3*a4_2)*idetM;
		Minv(0,1) = (-a1_2*a3_3*a4_4+a1_2*a3_4*a4_3+a1_3*a3_2*a4_4-a1_3*a3_4*a4_2-a1_4*a3_2*a4_3+a1_4*a3_3*a4_2)*idetM;
		Minv(0,2) = (a1_2*a2_3*a4_4-a1_2*a2_4*a4_3-a1_3*a2_2*a4_4+a1_3*a2_4*a4_2+a1_4*a2_2*a4_3-a1_4*a2_3*a4_2)*idetM;
		Minv(0,3) = (-a1_2*a2_3*a3_4+a1_2*a2_4*a3_3+a1_3*a2_2*a3_4-a1_3*a2_4*a3_2-a1_4*a2_2*a3_3+a1_4*a2_3*a3_2)*idetM;
		Minv(1,0) = (-a2_1*a3_3*a4_4+a2_1*a3_4*a4_3+a2_3*a3_1*a4_4-a2_3*a3_4*a4_1-a2_4*a3_1*a4_3+a2_4*a3_3*a4_1)*idetM;
		Minv(1,1) = (a1_1*a3_3*a4_4-a1_1*a3_4*a4_3-a1_3*a3_1*a4_4+a1_3*a3_4*a4_1+a1_4*a3_1*a4_3-a1_4*a3_3*a4_1)*idetM;
		Minv(1,2) = (-a1_1*a2_3*a4_4+a1_1*a2_4*a4_3+a1_3*a2_1*a4_4-a1_3*a2_4*a4_1-a1_4*a2_1*a4_3+a1_4*a2_3*a4_1)*idetM;
		Minv(1,3) = (a1_1*a2_3*a3_4-a1_1*a2_4*a3_3-a1_3*a2_1*a3_4+a1_3*a2_4*a3_1+a1_4*a2_1*a3_3-a1_4*a2_3*a3_1)*idetM;
		Minv(2,0) = (a2_1*a3_2*a4_4-a2_1*a3_4*a4_2-a2_2*a3_1*a4_4+a2_2*a3_4*a4_1+a2_4*a3_1*a4_2-a2_4*a3_2*a4_1)*idetM;
		Minv(2,1) = (-a1_1*a3_2*a4_4+a1_1*a3_4*a4_2+a1_2*a3_1*a4_4-a1_2*a3_4*a4_1-a1_4*a3_1*a4_2+a1_4*a3_2*a4_1)*idetM;
		Minv(2,2) = (a1_1*a2_2*a4_4-a1_1*a2_4*a4_2-a1_2*a2_1*a4_4+a1_2*a2_4*a4_1+a1_4*a2_1*a4_2-a1_4*a2_2*a4_1)*idetM;
		Minv(2,3) = (-a1_1*a2_2*a3_4+a1_1*a2_4*a3_2+a1_2*a2_1*a3_4-a1_2*a2_4*a3_1-a1_4*a2_1*a3_2+a1_4*a2_2*a3_1)*idetM;
		Minv(3,0) = (-a2_1*a3_2*a4_3+a2_1*a3_3*a4_2+a2_2*a3_1*a4_3-a2_2*a3_3*a4_1-a2_3*a3_1*a4_2+a2_3*a3_2*a4_1)*idetM;
		Minv(3,1) = (a1_1*a3_2*a4_3-a1_1*a3_3*a4_2-a1_2*a3_1*a4_3+a1_2*a3_3*a4_1+a1_3*a3_1*a4_2-a1_3*a3_2*a4_1)*idetM;
		Minv(3,2) = (-a1_1*a2_2*a4_3+a1_1*a2_3*a4_2+a1_2*a2_1*a4_3-a1_2*a2_3*a4_1-a1_3*a2_1*a4_2+a1_3*a2_2*a4_1)*idetM;
		Minv(3,3) = (a1_1*a2_2*a3_3-a1_1*a2_3*a3_2-a1_2*a2_1*a3_3+a1_2*a2_3*a3_1+a1_3*a2_1*a3_2-a1_3*a2_2*a3_1)*idetM;
	}
	else
		luma_implement;
  return Minv;
}

} // numerics

} // luma

#endif
