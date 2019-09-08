// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "LinearAlgebra/BlockLinAlg/MatrixBlock_2x2_InverseLU.h"
#include "LinearAlgebra/BlockLinAlg/VectorBlock_2.h"
#include "LinearAlgebra/DenseLinAlg/InverseLU.h"

#include "LinearAlgebra/DenseLinAlg/StaticSize/MatrixS.h"
#include "LinearAlgebra/DenseLinAlg/StaticSize/VectorS.h"
#include "LinearAlgebra/DenseLinAlg/DynamicSize/MatrixD.h"
#include "LinearAlgebra/DenseLinAlg/DynamicSize/VectorD.h"

//This computes computes a matrix inverse using LU decomposition without pivoting

namespace SANS
{
namespace BLA
{

template< class M00, class M01, class M10, class M11, class MatrixType >
void MatrixBlock_2x2_LUSolver<M00, M01, M10, M11, MatrixType>::
Solve( const FactorType& Factorized, const Real sgn, MatrixType& res )
{
  const MatrixBlock_2x2<M00, M01, M10, M11>& AFac = Factorized.MatrixFac;

  SANS_ASSERT( AFac.m00.m() == AFac.m00.n() );
  SANS_ASSERT( AFac.m11.m() == AFac.m11.n() );

  //Forward solve (L has non-one diagonals)
  M00 invdiag00(AFac.m00.size());
  invdiag00 = 0;
  invdiag00 = DLA::InverseLU::Inverse( AFac.m00 );
  res.scale_row0(invdiag00);
  res.axpy_rows01(-AFac.m10);

  M11 invdiag11(AFac.m11.size());
  invdiag11 = DLA::InverseLU::Inverse( AFac.m11 );
  res.scale_row1(invdiag11);

  //Backward solve (U has one diagonals, no scale needed)
  res.axpy_rows10(-AFac.m01);

  if ( sgn != 1 )
    res *= sgn;
}


//Explicit instantiations
typedef DLA::MatrixD<Real> MD_Real;

typedef DLA::MatrixD<DLA::MatrixS<2,2,Real>> MD_MS22;
typedef DLA::MatrixD<DLA::VectorS<2,Real>> MD_VecS2;

typedef DLA::MatrixD<DLA::MatrixD<Real>> MD_MD_Real;
typedef DLA::VectorD<DLA::VectorD<Real>> VecD_VecD_Real;

template struct MatrixBlock_2x2_LUSolver<MD_Real, MD_Real, MD_Real, MD_Real, VectorBlock_2<MD_Real, MD_Real>>;
template struct MatrixBlock_2x2_LUSolver<MD_Real, MD_Real, MD_Real, MD_Real, MatrixBlock_2x2<MD_Real, MD_Real, MD_Real, MD_Real>>;

template struct MatrixBlock_2x2_LUSolver<MD_MS22, MD_MS22, MD_MS22, MD_MS22, VectorBlock_2<MD_VecS2, MD_VecS2>>;
template struct MatrixBlock_2x2_LUSolver<MD_MS22, MD_MS22, MD_MS22, MD_MS22, VectorBlock_2<MD_MS22, MD_MS22>>;

template struct MatrixBlock_2x2_LUSolver<MD_MD_Real, MD_MD_Real, MD_MD_Real, MD_MD_Real, VectorBlock_2<VecD_VecD_Real, VecD_VecD_Real>>;

} //namespace BLA
} //namespace SANS
