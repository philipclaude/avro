// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MatMul_Native.h"

#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD.h"
#include "LinearAlgebra/DenseLinAlg/DynamicSize/VectorD.h"

#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"

namespace SANS
{
namespace DLA
{

template<class TL, class TR, class T>
void MatMul_Native<TL, TR, T>::value(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res )
{
  res = 0;
  plus(ML, MR, sgn, res);
}

template<class TL, class TR, class T>
void MatMul_Native<TL, TR, T>::plus(const MatrixDView<TL>& ML, const MatrixDView<TR>& MR, const Real sgn, MatrixDView<T>& res )
{
  //No need to do any optimization here as we are likely working with large sparse matrices anyways
  SANS_ASSERT( ML.m() == res.m() );
  SANS_ASSERT( MR.n() == res.n() );
  SANS_ASSERT( ML.n() == MR.m() );

  for (int i = 0; i < ML.m(); i++)
    for (int j = 0; j < ML.n(); j++)
      for (int k = 0; k < MR.n(); k++)
        res(i,k) += sgn*(ML(i,j)*MR(j,k));
}


template class MatMul_Native< MatrixD<Real>, MatrixD<Real>, MatrixD<Real> >;
template class MatMul_Native< MatrixD<Real>, VectorD<Real>, VectorD<Real> >;

template class MatMul_Native< MatrixD<MatrixS<1, 2, Real>>, VectorD<VectorS<2, Real>>, VectorD<Real> >;
template class MatMul_Native< MatrixD<MatrixS<1, 3, Real>>, VectorD<VectorS<3, Real>>, VectorD<Real> >;

template class MatMul_Native< MatrixD<MatrixS<2, 1, Real>>, VectorD<Real>, VectorD<VectorS<2, Real>> >;
template class MatMul_Native< MatrixD<MatrixS<2, 2, Real>>, VectorD<VectorS<2, Real>>, VectorD<VectorS<2, Real>> >;
template class MatMul_Native< MatrixD<MatrixS<2, 3, Real>>, VectorD<VectorS<3, Real>>, VectorD<VectorS<2, Real>> >;

template class MatMul_Native< MatrixD<MatrixS<3, 1, Real>>, VectorD<Real>, VectorD<VectorS<3, Real>> >;
template class MatMul_Native< MatrixD<MatrixS<3, 2, Real>>, VectorD<VectorS<2, Real>>, VectorD<VectorS<3, Real>> >;
template class MatMul_Native< MatrixD<MatrixS<3, 3, Real>>, VectorD<VectorS<3, Real>>, VectorD<VectorS<3, Real>> >;

template class MatMul_Native< MatrixD<MatrixS<1,3,Real>>, MatrixD<MatrixS<3,1,Real>>, MatrixD<Real> >;
template class MatMul_Native< MatrixD<Real>, VectorD<VectorS<3,Real> >, VectorD<VectorS<3,Real> > >;

template class MatMul_Native< MatrixS<1, 1, MatrixS<4, 4, Real> >,
                              VectorS<1, MatrixS<4, 4, Real> >,
                              MatrixS<4, 4, Real>>;

template class MatMul_Native< MatrixS<1, 1, MatrixS<4, 4, Real> >,
                              VectorS<1, Real>,
                              MatrixS<4, 4, Real>>;

template class MatMul_Native<VectorS<2, MatrixS<6, 6, Real> >, MatrixS<1, 1, VectorS<6, Real> >, VectorS<2, VectorS<6, Real> > >;
template class MatMul_Native<MatrixS<1, 2, MatrixS<6, 6, Real> >, VectorS<2, VectorS<6, Real> >, VectorS<6, Real> >;

template class MatMul_Native< MatrixS<3,3,Real>, MatrixS<3,1,Real>, MatrixS<3,1,Real> >;

template class MatMul_Native< Real, VectorS<10,Real>,VectorS<10,Real> >;

} //namespace DLA
} //namespace SANS
