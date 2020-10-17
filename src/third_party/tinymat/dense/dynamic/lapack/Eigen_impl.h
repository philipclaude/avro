// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(EIGEN_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include "Eigen.h"
#include "tools/SANSException.h"

//#ifdef DLA_LAPACK

#include <vector>

#include "tinymat/dense/dynamic/MatrixD.h"
#include "tinymat/dense/dynamic/MatrixD_Transpose.h"
#include "tinymat/dense/dynamic/VectorD.h"
#include "tinymat/dense/tools/dense_LAPACK.h"

//Eigen value/vector routines. The lapack routines for single/double are defined in
//as part of the explicit instantiation in the .cpp files

namespace tinymat
{
namespace DLA
{

// A is the matrix we want the eigen values from
// wr is the real part of the eigen values
// wi is the imaginary part of the eigen values
template<class T>
void LAPACK_Eigen<T>::Value( const MatrixDView<T>& A, VectorDView<T>& wr, VectorDView<T>& wi )
{
  SANS_ASSERT( A.m() == A.n()  );
  SANS_ASSERT( A.m() == wr.m() );
  SANS_ASSERT( A.m() == wi.m() );

//  lapack_int INFO = LAPACKE_GEEV( LAPACK_ROW_MAJOR, 'n', 'n',
//                                  A.m(), &A(0,0), A.stride(), &wr(0,0),
//                                  &wi(0,0), NULL, A.m(), NULL, A.m() );

  char jobvl = 'N'; // left eigenvector
  char jobvr = 'N'; // right eigenvector
  T *vl = NULL;     // left eigenvector place holder
  int ldvl = 1;
  T *vr = NULL;     // right eigenvector place holder
  int ldvr = 1;
  int m = A.m();
  int stride = A.stride();
  int lwork = 16*A.m();
  std::vector<T> work( lwork );
  int INFO;

  LAPACK_GEEV(&jobvl,&jobvr,&m,const_cast<T*>(&A(0,0)),&stride,&wr[0],&wi[0],vl,&ldvl,vr,&ldvr,&work[0],&lwork,&INFO);

  SANS_ASSERT_MSG( INFO == 0, "INFO == %d", INFO );
}

// A is the matrix we want the eigen vectors from
// vl is the left eigen vectors
// vr is the right eigen vectors
// The eigen vectors are returned such that
// A*vr == vr*diag(wr)
// vl*A == diag(wr)*vl
template<class T>
void LAPACK_Eigen<T>::Vectors( MatrixDView<T>& A, MatrixDView<T>& vl, MatrixDView<T>& vr )
{
  SANS_ASSERT( A.m() == A.n()  );
  SANS_ASSERT( vl.m() == vl.n() );
  SANS_ASSERT( vr.m() == vr.n() );
  SANS_ASSERT( A.m() == vl.m() );
  SANS_ASSERT( A.m() == vr.m() );

  VectorD<T> wr(A.m());
  VectorD<T> wi(A.m());

  MatrixD<T> vrt(vr.m(), vr.n());

//  lapack_int INFO = LAPACKE_GEEV( LAPACK_ROW_MAJOR, 'v', 'v',
//                                  A.m(), &A(0,0), A.stride(), &wr(0,0),
//                                  &wi(0,0), &vl(0,0), vl.m(), &vr(0,0), vr.m() );

  char jobvl = 'V'; // left eigenvector (which is r-eig of row-matrix)
  char jobvr = 'V'; // right eigenvector
  int ldvl = vl.stride();
  int ldvr = vrt.stride();
  int m = A.m();
  int stride = A.stride();
  int lwork = 16*m;
  std::vector<T> work( lwork );
  int INFO;

  // Note that left and right are "flipped" with respect to the Fortran lapack call
  // because we assume that A is row-major
  LAPACK_GEEV(&jobvl, &jobvr,
              &m, &A(0,0), &stride,
              &wr[0], &wi[0],
              &vrt(0,0), &ldvr,
              &vl(0,0), &ldvl,
              &work[0],&lwork,&INFO);

  vr = Transpose(vrt);

  SANS_ASSERT_MSG( INFO == 0, "INFO == %d", INFO );
}


// A is the matrix we want the eigen values/vectors from
// wr is the real part of the eigen values
// wi is the imaginary part of the eigen values
// vl is the left eigen vectors
// vr is the right eigen vectors
// The eigen vectors are returned such that
// A*vr == vr*diag(wr)
// vl*A == diag(wr)*vl
template<class T>
void LAPACK_Eigen<T>::System( MatrixDView<T>& A, VectorDView<T>& wr, VectorDView<T>& wi,
                                                 MatrixDView<T>& vl, MatrixDView<T>& vr )
{
  SANS_ASSERT( A.m() == A.n()  );
  SANS_ASSERT( vl.m() == vl.n() );
  SANS_ASSERT( vr.m() == vr.n() );
  SANS_ASSERT( A.m() == vl.m() );
  SANS_ASSERT( A.m() == vr.m() );

  MatrixD<T> tmp(A);
  MatrixD<T> vrt(vr.m(), vr.n());

//  lapack_int INFO = LAPACKE_GEEV( LAPACK_ROW_MAJOR, 'v', 'v',
//                                  A.m(), &A(0,0), A.stride(), &wr(0,0),
//                                  &wi(0,0), &vl(0,0), vl.m(), &vr(0,0), vr.m() );

  char jobvl = 'V'; // left eigenvector (which is r-eig of row-matrix)
  char jobvr = 'V'; // right eigenvector
  int ldvl = vl.stride();
  int ldvr = vrt.stride();
  int m = A.m();
  int stride = A.stride();
  int lwork = 16*A.m();
  std::vector<T> work( lwork );
  int INFO;

  // Note that left and right are "flipped" with respect to the Fortran lapack call
  // because we assume that A is row-major
  LAPACK_GEEV(&jobvl, &jobvr,
              &m, &tmp(0,0), &stride,
              &wr[0], &wi[0],
              &vrt(0,0), &ldvr,
              &vl(0,0), &ldvl,
              &work[0],&lwork,&INFO);

  vr = Transpose(vrt);

  SANS_ASSERT_MSG( INFO == 0, "INFO == %d", INFO );
}

} //namespace tinymat
} //namespace DLA

//#endif
