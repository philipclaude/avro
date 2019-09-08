// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MatrixBlock_2x2_Decompose_LU.h"
//#include "numpack/DenseLinAlg/DynamicSize/MatrixD.h"
#include "numpack/DenseLinAlg/StaticSize/MatrixS.h"

#include "numpack/DenseLinAlg/InverseLU.h"

//Perform a LU decomposition without pivoting and unit diagonal on U

namespace numpack 
{
namespace BLA
{

//-----------------------------------------------------------------------------
template< class M00, class M01, class M10, class M11 >
void MatrixBlock_2x2_LU<M00, M01, M10, M11>::Decompose( MatrixBlock_2x2<M00, M01, M10, M11>& Matrix )
{
  /*/Perform the LU decomposition
   * L = [M00,         0             ]     U = [I, inv(M00)*M01]
   *     [M10, M11 - M10*inv(M00)*M01]         [0,     I       ]
   */

  Matrix.m01 = DLA::InverseLU::Solve(Matrix.m00, Matrix.m01);
  Matrix.m11 -= Matrix.m10*Matrix.m01;
}


//Explicit instantiations
typedef DLA::MatrixD<Real> MD_Real;
typedef DLA::MatrixD<DLA::MatrixS<2,2,Real>> MD_MS2;

typedef DLA::MatrixD<DLA::MatrixD<Real>> MD_MD_Real;

template struct MatrixBlock_2x2_LU<MD_Real, MD_Real, MD_Real, MD_Real>;
template struct MatrixBlock_2x2_LU<MD_MS2, MD_MS2, MD_MS2, MD_MS2>;

template struct MatrixBlock_2x2_LU<MD_MD_Real, MD_MD_Real, MD_MD_Real, MD_MD_Real>;

} //namespace BLA
} //namespace numpack 
