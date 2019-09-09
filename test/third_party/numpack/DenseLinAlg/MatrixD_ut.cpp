// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "unit_tester.hpp"

#include "numpack/DenseLinAlg/DynamicSize/MatrixD.h"
#include "numpack/DenseLinAlg/DynamicSize/MatrixD_Sub.h"
#include "numpack/DenseLinAlg/StaticSize/MatrixS.h"

#include "chkMatrixD_btest.h"

#include <iostream>
using namespace numpack::DLA;

namespace numpack
{
namespace DLA
{
//Explicitly instantiate the class to generate all the functions so that coverage information is correct
template class MatrixD<Real>;

typedef MatrixDView<Real> Mat;
typedef OpAddD< Mat, Mat, true > MatAddT;

template class OpMulD< Mat, Mat >;
template class OpMulD< MatAddT, Mat >;
template class OpMulD< Mat, MatAddT >;
template class OpMulD< MatAddT, MatAddT >;

template class OpAddD< Mat, Mat, true >;
template class OpSubD< Mat, Mat, true >;
template class OpMulDScalar< Mat, true >;
template class OpMulDFactor< Mat, OpMulDScalar<Mat, true> >;
template class OpMulDFactor< OpMulDScalar<Mat, true>, Mat >;

template class OpAddD< Mat, Mat, false >;
template class OpSubD< Mat, Mat, false >;
template class OpMulDScalar< Mat, false >;
template class OpMulDFactor< Mat, OpMulDScalar<Mat, false> >;
template class OpMulDFactor< OpMulDScalar<Mat, false>, Mat >;

}
}

//############################################################################//
UT_TEST_SUITE( DenseLinAlg_MatrixD )

//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixD_ctor )
{

  MatrixD<Real> M1(2, 2);

  UT_ASSERT_EQUALS( 4, M1.size() );
  UT_ASSERT_EQUALS( 2, M1.m() );
  UT_ASSERT_EQUALS( 2, M1.n() );
  UT_ASSERT_EQUALS( 2, M1.stride() );

  M1(0,0) = 1;
  M1(0,1) = 0;
  M1(1,0) = 3;
  M1(1,1) = 2;

  UT_ASSERT_EQUALS( 1, M1(0,0) );
  UT_ASSERT_EQUALS( 0, M1(0,1) );
  UT_ASSERT_EQUALS( 3, M1(1,0) );
  UT_ASSERT_EQUALS( 2, M1(1,1) );

  MatrixD<Real> M2(2, 2);
  M2 = 0;
  M2(0,0) = 1;
  M2(0,1) = 42;
  M2(1,0) = 1929;
  M2(1,1) = 2;

  UT_ASSERT_EQUALS( 4, M2.size() );
  UT_ASSERT_EQUALS( 2, M2.m() );
  UT_ASSERT_EQUALS( 2, M2.n() );
  UT_ASSERT_EQUALS( 2, M2.stride() );

  UT_ASSERT_EQUALS(    1, M2(0,0) );
  UT_ASSERT_EQUALS(   42, M2(0,1) );
  UT_ASSERT_EQUALS( 1929, M2(1,0) );
  UT_ASSERT_EQUALS(    2, M2(1,1) );

  Real v3[4] = {1, 42, 1929, 2};
  MatrixDView<Real> M3(v3, 2, 2);
  MatrixD<Real> M4( M3 );

  UT_ASSERT_EQUALS( 4, M4.size() );
  UT_ASSERT_EQUALS( 2, M4.m() );
  UT_ASSERT_EQUALS( 2, M4.n() );
  UT_ASSERT_EQUALS( 2, M4.stride() );

  UT_ASSERT_EQUALS(    1, M4(0,0) );
  UT_ASSERT_EQUALS(   42, M4(0,1) );
  UT_ASSERT_EQUALS( 1929, M4(1,0) );
  UT_ASSERT_EQUALS(    2, M4(1,1) );

  MatrixD<Real> M5( M4 );

  UT_ASSERT_EQUALS( 4, M5.size() );
  UT_ASSERT_EQUALS( 2, M5.m() );
  UT_ASSERT_EQUALS( 2, M5.n() );
  UT_ASSERT_EQUALS( 2, M5.stride() );

  UT_ASSERT_EQUALS(    1, M5(0,0) );
  UT_ASSERT_EQUALS(   42, M5(0,1) );
  UT_ASSERT_EQUALS( 1929, M5(1,0) );
  UT_ASSERT_EQUALS(    2, M5(1,1) );

  MatrixD<Real> M6 = M5;

  UT_ASSERT_EQUALS( 4, M6.size() );
  UT_ASSERT_EQUALS( 2, M6.m() );
  UT_ASSERT_EQUALS( 2, M6.n() );
  UT_ASSERT_EQUALS( 2, M6.stride() );

  UT_ASSERT_EQUALS(    1, M6(0,0) );
  UT_ASSERT_EQUALS(   42, M6(0,1) );
  UT_ASSERT_EQUALS( 1929, M6(1,0) );
  UT_ASSERT_EQUALS(    2, M6(1,1) );

  MatrixD<Real> M7 = M3;

  UT_ASSERT_EQUALS( 4, M7.size() );
  UT_ASSERT_EQUALS( 2, M7.m() );
  UT_ASSERT_EQUALS( 2, M7.n() );
  UT_ASSERT_EQUALS( 2, M7.stride() );

  UT_ASSERT_EQUALS(    1, M7(0,0) );
  UT_ASSERT_EQUALS(   42, M7(0,1) );
  UT_ASSERT_EQUALS( 1929, M7(1,0) );
  UT_ASSERT_EQUALS(    2, M7(1,1) );

  M7 = Identity();

  UT_ASSERT_EQUALS(  1, M7(0,0) );
  UT_ASSERT_EQUALS(  0, M7(0,1) );
  UT_ASSERT_EQUALS(  0, M7(1,0) );
  UT_ASSERT_EQUALS(  1, M7(1,1) );

  MatrixD<Real> M8(2,3);
  M8 = Identity();

  UT_ASSERT_EQUALS( 1, M8(0,0) );
  UT_ASSERT_EQUALS( 0, M8(0,1) );
  UT_ASSERT_EQUALS( 0, M8(0,2) );
  UT_ASSERT_EQUALS( 0, M8(1,0) );
  UT_ASSERT_EQUALS( 1, M8(1,1) );
  UT_ASSERT_EQUALS( 0, M8(1,2) );

  MatrixD<Real> M9({{1929,42},{-1,-3}});

  UT_ASSERT_EQUALS( 1929, M9(0,0) );
  UT_ASSERT_EQUALS(   42, M9(0,1) );
  UT_ASSERT_EQUALS(   -1, M9(1,0) );
  UT_ASSERT_EQUALS(   -3, M9(1,1) );

  MatrixD<Real> M10(2,2,0.);

  UT_ASSERT_EQUALS( 0, M10(0,0) );
  UT_ASSERT_EQUALS( 0, M10(0,1) );
  UT_ASSERT_EQUALS( 0, M10(1,0) );
  UT_ASSERT_EQUALS( 0, M10(1,1) );

  MatrixD< MatrixS<1,1,Real> > M11(2,2,0.);

  UT_ASSERT_EQUALS( 0, M11(0,0)(0,0) );
  UT_ASSERT_EQUALS( 0, M11(0,1)(0,0) );
  UT_ASSERT_EQUALS( 0, M11(1,0)(0,0) );
  UT_ASSERT_EQUALS( 0, M11(1,1)(0,0) );

  MatrixD< MatrixD<Real> > M12 = {{{2,2},{2,3}},
                                  {{3,2},{3,3}}};

  UT_ASSERT_EQUALS( 2, M12(0,0).m() );
  UT_ASSERT_EQUALS( 2, M12(0,0).n() );

  UT_ASSERT_EQUALS( 2, M12(0,1).m() );
  UT_ASSERT_EQUALS( 3, M12(0,1).n() );

  UT_ASSERT_EQUALS( 3, M12(1,0).m() );
  UT_ASSERT_EQUALS( 2, M12(1,0).n() );

  UT_ASSERT_EQUALS( 3, M12(1,1).m() );
  UT_ASSERT_EQUALS( 3, M12(1,1).n() );

}
UT_TEST_CASE_END( MatrixD_ctor)

//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixD_Assign )
{
  MatrixD<Real> M1(2, 2);
  M1 = 0;
  M1(0,0) = 1;
  M1(1,0) = 3;
  M1(1,1) = 2;

  MatrixD<Real> M2(2, 2);

  M2 = M1;

  UT_ASSERT_EQUALS( 1, M2(0,0) );
  UT_ASSERT_EQUALS( 0, M2(0,1) );
  UT_ASSERT_EQUALS( 3, M2(1,0) );
  UT_ASSERT_EQUALS( 2, M2(1,1) );

  M2 = -M1;

  UT_ASSERT_EQUALS( -1, M2(0,0) );
  UT_ASSERT_EQUALS(  0, M2(0,1) );
  UT_ASSERT_EQUALS( -3, M2(1,0) );
  UT_ASSERT_EQUALS( -2, M2(1,1) );

  //Make sure that an exception is thrown when dimensions do not match
  MatrixD<Real> M3(2, 3);

  //UT_ASSERT_THROW( M3 = M1, AssertionException );

  MatrixD<Real> M4(3, 2);

  //UT_ASSERT_THROW( M4 = M1, AssertionException );
}
UT_TEST_CASE_END( MatrixD_Assign )


//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixD_CompoundAssign )
{
  MatrixD<Real> M1(2, 2), M2(2, 2);
  M1 = 0;
  M1(0,0) = 1;
  M1(1,0) = 3;
  M1(1,1) = 2;

  M2 = 0;
  M2(0,0) = 3;
  M2(0,1) = 1;
  M2(1,0) = 2;
  M2(1,1) = 4;

  M2 += M1;

  UT_ASSERT_EQUALS( 1, M1(0,0) );
  UT_ASSERT_EQUALS( 0, M1(0,1) );
  UT_ASSERT_EQUALS( 3, M1(1,0) );
  UT_ASSERT_EQUALS( 2, M1(1,1) );

  UT_ASSERT_EQUALS( 4, M2(0,0) );
  UT_ASSERT_EQUALS( 1, M2(0,1) );
  UT_ASSERT_EQUALS( 5, M2(1,0) );
  UT_ASSERT_EQUALS( 6, M2(1,1) );

  M1 -= M2;

  UT_ASSERT_EQUALS( -3, M1(0,0) );
  UT_ASSERT_EQUALS( -1, M1(0,1) );
  UT_ASSERT_EQUALS( -2, M1(1,0) );
  UT_ASSERT_EQUALS( -4, M1(1,1) );

  UT_ASSERT_EQUALS( 4, M2(0,0) );
  UT_ASSERT_EQUALS( 1, M2(0,1) );
  UT_ASSERT_EQUALS( 5, M2(1,0) );
  UT_ASSERT_EQUALS( 6, M2(1,1) );

  M2 *= 2;

  UT_ASSERT_EQUALS(  8, M2(0,0) );
  UT_ASSERT_EQUALS(  2, M2(0,1) );
  UT_ASSERT_EQUALS( 10, M2(1,0) );
  UT_ASSERT_EQUALS( 12, M2(1,1) );

  M2 /= 2;

  UT_ASSERT_EQUALS( 4, M2(0,0) );
  UT_ASSERT_EQUALS( 1, M2(0,1) );
  UT_ASSERT_EQUALS( 5, M2(1,0) );
  UT_ASSERT_EQUALS( 6, M2(1,1) );

  +M2; //This should do nothing

  UT_ASSERT_EQUALS( 4, M2(0,0) );
  UT_ASSERT_EQUALS( 1, M2(0,1) );
  UT_ASSERT_EQUALS( 5, M2(1,0) );
  UT_ASSERT_EQUALS( 6, M2(1,1) );


  //Make sure that an exception is thrown when dimensions do not match
  MatrixD<Real> M3(2, 3);
  MatrixD<Real> M4(3, 2);

  //UT_ASSERT_THROW( M3 += M1, AssertionException );
  //UT_ASSERT_THROW( M4 += M1, AssertionException );

  //UT_ASSERT_THROW( M3 -= M1, AssertionException );
  //UT_ASSERT_THROW( M4 -= M1, AssertionException );

}
UT_TEST_CASE_END( MatrixD_CompoundAssign )


//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixD_Add )
{
  MatrixD<Real> M1(2, 2), M2(2, 2), M3(2, 2);
  M1 = 0;
  M1(0,0) = 1;
  M1(1,0) = 3;
  M1(1,1) = 2;

  M2 = 0;
  M2(0,0) = 3;
  M2(0,1) = 1;
  M2(1,0) = 2;
  M2(1,1) = 4;

  M3 = M2 + M1;

  UT_ASSERT_EQUALS( 1, M1(0,0) );
  UT_ASSERT_EQUALS( 0, M1(0,1) );
  UT_ASSERT_EQUALS( 3, M1(1,0) );
  UT_ASSERT_EQUALS( 2, M1(1,1) );

  UT_ASSERT_EQUALS( 3, M2(0,0) );
  UT_ASSERT_EQUALS( 1, M2(0,1) );
  UT_ASSERT_EQUALS( 2, M2(1,0) );
  UT_ASSERT_EQUALS( 4, M2(1,1) );

  UT_ASSERT_EQUALS( 4, M3(0,0) );
  UT_ASSERT_EQUALS( 1, M3(0,1) );
  UT_ASSERT_EQUALS( 5, M3(1,0) );
  UT_ASSERT_EQUALS( 6, M3(1,1) );

  M3 += M2 + M1;

  UT_ASSERT_EQUALS( 2*4, M3(0,0) );
  UT_ASSERT_EQUALS( 2*1, M3(0,1) );
  UT_ASSERT_EQUALS( 2*5, M3(1,0) );
  UT_ASSERT_EQUALS( 2*6, M3(1,1) );

  M3 = M1 - M2;

  UT_ASSERT_EQUALS( -2, M3(0,0) );
  UT_ASSERT_EQUALS( -1, M3(0,1) );
  UT_ASSERT_EQUALS(  1, M3(1,0) );
  UT_ASSERT_EQUALS( -2, M3(1,1) );

  M3 += M1 - M2;

  UT_ASSERT_EQUALS( -2*2, M3(0,0) );
  UT_ASSERT_EQUALS( -1*2, M3(0,1) );
  UT_ASSERT_EQUALS(  1*2, M3(1,0) );
  UT_ASSERT_EQUALS( -2*2, M3(1,1) );

  M3 -= M1 - M2;

  UT_ASSERT_EQUALS( -2, M3(0,0) );
  UT_ASSERT_EQUALS( -1, M3(0,1) );
  UT_ASSERT_EQUALS(  1, M3(1,0) );
  UT_ASSERT_EQUALS( -2, M3(1,1) );

  M1(0,0) = 1;
  M1(0,1) = 0;
  M1(1,0) = 3;
  M1(1,1) = 2;

  M2(0,0) = 3;
  M2(0,1) = 1;
  M2(1,0) = 2;
  M2(1,1) = 4;

  M3(0,0) = 4;
  M3(0,1) = 3;
  M3(1,0) = 2;
  M3(1,1) = 1;

  MatrixD<Real> M4(2, 2);

  M4 = M3 + M2 + M1;

  UT_ASSERT_EQUALS( 8, M4(0,0) );
  UT_ASSERT_EQUALS( 4, M4(0,1) );
  UT_ASSERT_EQUALS( 7, M4(1,0) );
  UT_ASSERT_EQUALS( 7, M4(1,1) );

  M4 = M3 - M2 + M1;

  UT_ASSERT_EQUALS(  2, M4(0,0) );
  UT_ASSERT_EQUALS(  2, M4(0,1) );
  UT_ASSERT_EQUALS(  3, M4(1,0) );
  UT_ASSERT_EQUALS( -1, M4(1,1) );

  M4 = M3 + M2 - M1;

  UT_ASSERT_EQUALS(  6, M4(0,0) );
  UT_ASSERT_EQUALS(  4, M4(0,1) );
  UT_ASSERT_EQUALS(  1, M4(1,0) );
  UT_ASSERT_EQUALS(  3, M4(1,1) );

  //It is ok for what is on the left to be on the right as long as it is the first thing
  M4 = M4 + M3 + M2 - M1;
  UT_ASSERT_EQUALS(  2*6, M4(0,0) );
  UT_ASSERT_EQUALS(  2*4, M4(0,1) );
  UT_ASSERT_EQUALS(  2*1, M4(1,0) );
  UT_ASSERT_EQUALS(  2*3, M4(1,1) );

  //The same object can be on both left and right as long is only addition/subrtraction and scalar multiplication is involved
  M4 = M3 + M4;
  UT_ASSERT_EQUALS(  2*6+4, M4(0,0) );
  UT_ASSERT_EQUALS(  2*4+3, M4(0,1) );
  UT_ASSERT_EQUALS(  2*1+2, M4(1,0) );
  UT_ASSERT_EQUALS(  2*3+1, M4(1,1) );

  M4 = M3 - M4;
  UT_ASSERT_EQUALS( -2*6, M4(0,0) );
  UT_ASSERT_EQUALS( -2*4, M4(0,1) );
  UT_ASSERT_EQUALS( -2*1, M4(1,0) );
  UT_ASSERT_EQUALS( -2*3, M4(1,1) );

  //Make sure that an exception is thrown when dimensions do not match
  MatrixD<Real> M23(2, 3);
  MatrixD<Real> M32(3, 2);

  //UT_ASSERT_THROW( M23 =  M2 +  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 = M23 +  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 =  M1 + M23, AssertionException );

  //UT_ASSERT_THROW( M32 =  M2 +  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 = M32 +  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 =  M1 + M32, AssertionException );

  //UT_ASSERT_THROW( M23 =  M2 -  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 = M23 -  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 =  M1 - M23, AssertionException );

  //UT_ASSERT_THROW( M32 =  M2 -  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 = M32 -  M1, AssertionException );
  //UT_ASSERT_THROW(  M3 =  M1 - M32, AssertionException );

}
UT_TEST_CASE_END( MatrixD_Add )


//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixD_ops1 )
{
  const Real data = 3;
  MatrixD<Real> m1(1, 1);
  MatrixD<Real> m2(1, 1);
  MatrixD<Real> m3(1, 1), m4(1, 1), m5(1, 1);

  m1 = data;
  m2 = data;

  // ctors
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );

  // assignment
  m3 = m1;
  m4 = data;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m3(0,0) );
  UT_ASSERT_EQUALS(  3, m4(0,0) );

  m1 = m2 = m3 = data;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  3, m3(0,0) );

  m4 = data;
  m1 = m2 = m3 = m4;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  3, m3(0,0) );
  UT_ASSERT_EQUALS(  3, m4(0,0) );

  // unary
  m2 = +m1;
  m3 = -m1;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS( -3, m3(0,0) );

  // binary accumulation
  m3 = m1;
  m3 *= data;
  UT_ASSERT_EQUALS(  9, m3(0,0) );

  m3 = data;
  m3 *= m1;
  UT_ASSERT_EQUALS(  9, m3(0,0) );

  m1 = data;
  m2 = m1;
  m3 = m1;
  m2 += m1;
  m3 -= m1;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  6, m2(0,0) );
  UT_ASSERT_EQUALS(  0, m3(0,0) );

  // binary operators
  m1 = data;
  //m2 = m1 + data;
  //m3 = m1 - data;
  m4 = m1 * data;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  //UT_ASSERT_EQUALS(  6, m2(0,0) );
  //UT_ASSERT_EQUALS(  0, m3(0,0) );
  UT_ASSERT_EQUALS(  9, m4(0,0) );

  //m2 = data + m1;
  //m3 = data - m1;
  m4 = data * m1;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  //UT_ASSERT_EQUALS(  6, m2(0,0) );
  //UT_ASSERT_EQUALS(  0, m3(0,0) );
  UT_ASSERT_EQUALS(  9, m4(0,0) );

  m1 = m2 = data;
  m3 = m1 + m2;
  m4 = m1 - m2;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  6, m3(0,0) );
  UT_ASSERT_EQUALS(  0, m4(0,0) );

  // arithmetic combinations

  m1 = m2 = data;
  m3 = m1 + m2;
  m4 = m1 + m2 + m3;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  6, m3(0,0) );
  UT_ASSERT_EQUALS( 12, m4(0,0) );

  m2 += m1;
  m3 += m1 + m2;
  m4 += m1 + m2 + m3;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  6, m2(0,0) );
  UT_ASSERT_EQUALS( 15, m3(0,0) );
  UT_ASSERT_EQUALS( 36, m4(0,0) );

  m3 = m1 - m2;
  m4 = m1 - m2 - m3;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  6, m2(0,0) );
  UT_ASSERT_EQUALS( -3, m3(0,0) );
  UT_ASSERT_EQUALS(  0, m4(0,0) );

  m2 -= m1;
  m3 -= m1 - m2;
  m4 -= m1 - m2 - m3;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS( -3, m3(0,0) );
  UT_ASSERT_EQUALS( -3, m4(0,0) );

  m3 = m1 - m2;
  m4 = m1 + m2 - m3;
  m5 = m1 - m2 + m3;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  0, m3(0,0) );
  UT_ASSERT_EQUALS(  6, m4(0,0) );
  UT_ASSERT_EQUALS(  0, m5(0,0) );

  m5 = (m1 + m2) + (m3 + m4);
  UT_ASSERT_EQUALS( 12, m5(0,0) );
  m5 = (m1 + m2) + (m3 - m4);
  UT_ASSERT_EQUALS(  0, m5(0,0) );
  m5 = (m1 + m2) - (m3 + m4);
  UT_ASSERT_EQUALS(  0, m5(0,0) );
  m5 = (m1 + m2) - (m3 - m4);
  UT_ASSERT_EQUALS( 12, m5(0,0) );
  m5 = (m1 - m2) + (m3 + m4);
  UT_ASSERT_EQUALS(  6, m5(0,0) );
  m5 = (m1 - m2) + (m3 - m4);
  UT_ASSERT_EQUALS( -6, m5(0,0) );
  m5 = (m1 - m2) - (m3 + m4);
  UT_ASSERT_EQUALS( -6, m5(0,0) );
  m5 = (m1 - m2) - (m3 - m4);
  UT_ASSERT_EQUALS(  6, m5(0,0) );
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  0, m3(0,0) );
  UT_ASSERT_EQUALS(  6, m4(0,0) );

  m5 += (m1 + m2) + (m3 + m4);
  m5 += (m1 + m2) + (m3 - m4);
  m5 += (m1 + m2) - (m3 + m4);
  m5 += (m1 + m2) - (m3 - m4);
  m5 += (m1 - m2) + (m3 + m4);
  m5 += (m1 - m2) + (m3 - m4);
  m5 += (m1 - m2) - (m3 + m4);
  m5 += (m1 - m2) - (m3 - m4);
  UT_ASSERT_EQUALS( 30, m5(0,0) );
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  0, m3(0,0) );
  UT_ASSERT_EQUALS(  6, m4(0,0) );

  m5 -= (m1 + m2) + (m3 + m4);
  m5 -= (m1 + m2) + (m3 - m4);
  m5 -= (m1 + m2) - (m3 + m4);
  m5 -= (m1 + m2) - (m3 - m4);
  m5 -= (m1 - m2) + (m3 + m4);
  m5 -= (m1 - m2) + (m3 - m4);
  m5 -= (m1 - m2) - (m3 + m4);
  m5 -= (m1 - m2) - (m3 - m4);
  UT_ASSERT_EQUALS(  6, m5(0,0) );
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  0, m3(0,0) );
  UT_ASSERT_EQUALS(  6, m4(0,0) );

  m1 = data;

  m2 = 4*m1;
  m3 = m2*7;
  UT_ASSERT_EQUALS(   3, m1(0,0) );
  UT_ASSERT_EQUALS(  12, m2(0,0) );
  UT_ASSERT_EQUALS(  84, m3(0,0) );

  m2 = m3/7;
  m1 = m2/4;
  UT_ASSERT_EQUALS(   3, m1(0,0) );
  UT_ASSERT_EQUALS(  12, m2(0,0) );
  UT_ASSERT_EQUALS(  84, m3(0,0) );

  m2 += 4*m1;
  m3 += m2*7;
  UT_ASSERT_EQUALS(   3, m1(0,0) );
  UT_ASSERT_EQUALS(  24, m2(0,0) );
  UT_ASSERT_EQUALS( 252, m3(0,0) );

  m2 -= 4*m1;
  m3 -= m2*7;
  UT_ASSERT_EQUALS(   3, m1(0,0) );
  UT_ASSERT_EQUALS(  12, m2(0,0) );
  UT_ASSERT_EQUALS( 168, m3(0,0) );

  m2 = 2*m1/7;
  m3 = 2*(m2*4);
  UT_ASSERT_EQUALS(             3, m1(0,0) );
  UT_ASSERT_EQUALS(      2.*3./7., m2(0,0) );
  UT_ASSERT_EQUALS( 2*2.*3./7.*4., m3(0,0) );

  //Check multiplication with 1x1 matricies
  m1 = 2;
  m2 = 3;
  m3 = m1*m2;
  UT_ASSERT_EQUALS( 6, m3(0,0) );
  m3 += m1*m2;
  UT_ASSERT_EQUALS( 2*6, m3(0,0) );
  m3 -= m1*m2;
  UT_ASSERT_EQUALS( 6, m3(0,0) );

  m1 = data;
  m2 = 4*data;
  m3 = 168;

  m5 = 2*(m1 + m2) + (m3 + m4)*3;
  UT_ASSERT_EQUALS(  552, m5(0,0) );
  m5 = 2*(m1 + m2) + (m3 - m4)*3;
  UT_ASSERT_EQUALS(  516, m5(0,0) );
  m5 = 2*(m1 + m2) - (m3 + m4)*3;
  UT_ASSERT_EQUALS( -492, m5(0,0) );
  m5 = 2*(m1 + m2) - (m3 - m4)*3;
  UT_ASSERT_EQUALS( -456, m5(0,0) );
  m5 = 2*(m1 - m2) + (m3 + m4)*3;
  UT_ASSERT_EQUALS(  504, m5(0,0) );
  m5 = 2*(m1 - m2) + (m3 - m4)*3;
  UT_ASSERT_EQUALS(  468, m5(0,0) );
  m5 = 2*(m1 - m2) - (m3 + m4)*3;
  UT_ASSERT_EQUALS( -540, m5(0,0) );
  m5 = 2*(m1 - m2) - (m3 - m4)*3;
  UT_ASSERT_EQUALS( -504, m5(0,0) );
  UT_ASSERT_EQUALS(    3, m1(0,0) );
  UT_ASSERT_EQUALS(   12, m2(0,0) );
  UT_ASSERT_EQUALS(  168, m3(0,0) );
  UT_ASSERT_EQUALS(    6, m4(0,0) );

  m5 += 2*(m1 + m2) + (m3 + m4)*3;
  UT_ASSERT_EQUALS(   48, m5(0,0) );
  m5 += 2*(m1 + m2) + (m3 - m4)*3;
  UT_ASSERT_EQUALS(  564, m5(0,0) );
  m5 += 2*(m1 + m2) - (m3 + m4)*3;
  UT_ASSERT_EQUALS(   72, m5(0,0) );
  m5 += 2*(m1 + m2) - (m3 - m4)*3;
  UT_ASSERT_EQUALS( -384, m5(0,0) );
  m5 += 2*(m1 - m2) + (m3 + m4)*3;
  UT_ASSERT_EQUALS(  120, m5(0,0) );
  m5 += 2*(m1 - m2) + (m3 - m4)*3;
  UT_ASSERT_EQUALS(  588, m5(0,0) );
  m5 += 2*(m1 - m2) - (m3 + m4)*3;
  UT_ASSERT_EQUALS(   48, m5(0,0) );
  m5 += 2*(m1 - m2) - (m3 - m4)*3;
  UT_ASSERT_EQUALS( -456, m5(0,0) );
  UT_ASSERT_EQUALS(    3, m1(0,0) );
  UT_ASSERT_EQUALS(   12, m2(0,0) );
  UT_ASSERT_EQUALS(  168, m3(0,0) );
  UT_ASSERT_EQUALS(    6, m4(0,0) );

  m5 -= 2*(m1 + m2) + (m3 + m4)*3;
  UT_ASSERT_EQUALS( -1008, m5(0,0) );
  m5 -= 2*(m1 + m2) + (m3 - m4)*3;
  UT_ASSERT_EQUALS( -1524, m5(0,0) );
  m5 -= 2*(m1 + m2) - (m3 + m4)*3;
  UT_ASSERT_EQUALS( -1032, m5(0,0) );
  m5 -= 2*(m1 + m2) - (m3 - m4)*3;
  UT_ASSERT_EQUALS(  -576, m5(0,0) );
  m5 -= 2*(m1 - m2) + (m3 + m4)*3;
  UT_ASSERT_EQUALS( -1080, m5(0,0) );
  m5 -= 2*(m1 - m2) + (m3 - m4)*3;
  UT_ASSERT_EQUALS( -1548, m5(0,0) );
  m5 -= 2*(m1 - m2) - (m3 + m4)*3;
  UT_ASSERT_EQUALS( -1008, m5(0,0) );
  m5 -= 2*(m1 - m2) - (m3 - m4)*3;
  UT_ASSERT_EQUALS(  -504, m5(0,0) );
  UT_ASSERT_EQUALS(     3, m1(0,0) );
  UT_ASSERT_EQUALS(    12, m2(0,0) );
  UT_ASSERT_EQUALS(   168, m3(0,0) );
  UT_ASSERT_EQUALS(     6, m4(0,0) );

  m5 = 2*(m1 + m2)*3;
  UT_ASSERT_EQUALS( 90, m5(0,0) );
  m5 = 2*3*(m1 + m2);
  UT_ASSERT_EQUALS( 90, m5(0,0) );
  m5 = (m1 + m2)*2*3;
  UT_ASSERT_EQUALS( 90, m5(0,0) );
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS( 12, m2(0,0) );

  m2 = +m1;
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  m3 = -m2;
  UT_ASSERT_EQUALS( -3, m3(0,0) );
  m4 = +(m1 + m2);
  UT_ASSERT_EQUALS(  6, m4(0,0) );
  m4 = +(m1 - m2);
  UT_ASSERT_EQUALS(  0, m4(0,0) );
  m4 = -(m1 + m2);
  UT_ASSERT_EQUALS( -6, m4(0,0) );
  m4 = -(m1 - m2);
  UT_ASSERT_EQUALS(  0, m4(0,0) );
  m4 = +(m1 + m2) + m3;
  UT_ASSERT_EQUALS(  3, m4(0,0) );
  m4 = -(m1 + m2) + m3;
  UT_ASSERT_EQUALS( -9, m4(0,0) );
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS( -3, m3(0,0) );

  m4 = +5*m1;
  UT_ASSERT_EQUALS(  15, m4(0,0) );
  m4 = -5*m1;
  UT_ASSERT_EQUALS( -15, m4(0,0) );
  m4 = +m1*5;
  UT_ASSERT_EQUALS(  15, m4(0,0) );
  m4 = -m1*5;
  UT_ASSERT_EQUALS( -15, m4(0,0) );
  m4 = +(5*m1);
  UT_ASSERT_EQUALS(  15, m4(0,0) );
  m4 = -(5*m1);
  UT_ASSERT_EQUALS( -15, m4(0,0) );
  m4 = +(m1*5);
  UT_ASSERT_EQUALS(  15, m4(0,0) );
  m4 = -(m1*5);
  UT_ASSERT_EQUALS( -15, m4(0,0) );
  UT_ASSERT_EQUALS(   3, m1(0,0) );

  int i0 = m1*m1;
  UT_ASSERT_EQUALS( 9, i0 );

  int i1 = (m1+m1)*m1;
  UT_ASSERT_EQUALS( 18, i1 );

  int i2 = m1*(m1+m1);
  UT_ASSERT_EQUALS( 18, i2 );

  int i3 = (m1+m1)*(m1+m1);
  UT_ASSERT_EQUALS( 36, i3 );
}
UT_TEST_CASE_END( MatrixD_ops1 )


//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixD_ops2 )
{

  MatrixD<Real> m1(2, 2);
  MatrixD<Real> m2(2, 2);
  MatrixD<Real> m3(2, 2), m4(2, 2), m5(2, 2);

  m1(0,0) = 1; m1(0,1) = 2;
  m1(1,0) = 3; m1(1,1) = 4;

  m2 = m1;

  // ctors
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );

  // assignment
  m3 = m1;
  m4 = 5;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m4, 5,5,5,5 ) );

  m2 = m3 = 3;
  UT_ASSERT( chkMatrixD22( m2, 3,3,3,3 ) );
  UT_ASSERT( chkMatrixD22( m3, 3,3,3,3 ) );

//  m3 = m2 = {1, 2,
//             3, 4};
//  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
//  UT_ASSERT( chkMatrixD22( m3, 1,2,3,4 ) );

  // unary
  m2 = +m1;
  m3 = -m1;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, -1,-2,-3,-4 ) );

  // binary accumulation
  m3 = m1;
  m4 = m1;
  m4 *= 5;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m4, 5,10,15,20 ) );

  m2 = 5;
  m3 = 5;
  m2 += m1;
  m3 -= m1;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 6,7,8,9 ) );
  UT_ASSERT( chkMatrixD22( m3, 4,3,2,1 ) );

  // binary operators
  //m2 = m1 + 3;
  //m3 = m1 - 3;
  m4 = m1 * 3;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  //UT_ASSERT( chkMatrixD22( m2, 4,5,6,7 ) );
  //UT_ASSERT( chkMatrixD22( m3, -2,-1,0,1 ) );
  UT_ASSERT( chkMatrixD22( m4, 3,6,9,12 ) );

  //m2 = 3 + m1;
  //m3 = 3 - m1;
  m4 = 3 * m1;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  //UT_ASSERT( chkMatrixD22( m2, 4,5,6,7 ) );
  //UT_ASSERT( chkMatrixD22( m3, 2,1,0,-1 ) );
  UT_ASSERT( chkMatrixD22( m4, 3,6,9,12 ) );

  m2 = 3;
  m3 = m1 + m2;
  m4 = m1 - m2;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,3,3,3 ) );
  UT_ASSERT( chkMatrixD22( m3, 4,5,6,7 ) );
  UT_ASSERT( chkMatrixD22( m4, -2,-1,0,1 ) );

  // arithmetic combinations

  m2 = m1;
  m3 = m1 + m2;
  m4 = m1 + m2 + m3;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 2,4,6,8 ) );
  UT_ASSERT( chkMatrixD22( m4, 4,8,12,16 ) );

  m2 += m1;
  m3 += m1 + m2;
  m4 += m1 + m2 + m3;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 2,4,6,8 ) );
  UT_ASSERT( chkMatrixD22( m3, 5,10,15,20 ) );
  UT_ASSERT( chkMatrixD22( m4, 12,24,36,48 ) );

  m3 = m1 - m2;
  m4 = m1 - m2 - m3;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 2,4,6,8 ) );
  UT_ASSERT( chkMatrixD22( m3, -1,-2,-3,-4 ) );
  UT_ASSERT( chkMatrixD22( m4, 0,0,0,0 ) );

  m2 -= m1;
  m3 -= m1 - m2;
  m4 -= m1 - m2 - m3;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, -1,-2,-3,-4 ) );
  UT_ASSERT( chkMatrixD22( m4, -1,-2,-3,-4 ) );

  m3 = m1 - m2;
  m4 = m1 + m2 - m3;
  m5 = m1 - m2 + m3;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );
  UT_ASSERT( chkMatrixD22( m5, 0,0,0,0 ) );

  m5 = (m1 + m2) + (m3 + m4);
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  m5 = (m1 + m2) + (m3 - m4);
  UT_ASSERT( chkMatrixD22( m5, 0,0,0,0 ) );
  m5 = (m1 + m2) - (m3 + m4);
  UT_ASSERT( chkMatrixD22( m5, 0,0,0,0 ) );
  m5 = (m1 + m2) - (m3 - m4);
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  m5 = (m1 - m2) + (m3 + m4);
  UT_ASSERT( chkMatrixD22( m5, 2,4,6,8 ) );
  m5 = (m1 - m2) + (m3 - m4);
  UT_ASSERT( chkMatrixD22( m5, -2,-4,-6,-8 ) );
  m5 = (m1 - m2) - (m3 + m4);
  UT_ASSERT( chkMatrixD22( m5, -2,-4,-6,-8 ) );
  m5 = (m1 - m2) - (m3 - m4);
  UT_ASSERT( chkMatrixD22( m5, 2,4,6,8 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );

  m5 += (m1 + m2) + (m3 + m4);
  m5 += (m1 + m2) + (m3 - m4);
  m5 += (m1 + m2) - (m3 + m4);
  m5 += (m1 + m2) - (m3 - m4);
  m5 += (m1 - m2) + (m3 + m4);
  m5 += (m1 - m2) + (m3 - m4);
  m5 += (m1 - m2) - (m3 + m4);
  m5 += (m1 - m2) - (m3 - m4);
  UT_ASSERT( chkMatrixD22( m5, 10,20,30,40 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );

  m5 -= (m1 + m2) + (m3 + m4);
  m5 -= (m1 + m2) + (m3 - m4);
  m5 -= (m1 + m2) - (m3 + m4);
  m5 -= (m1 + m2) - (m3 - m4);
  m5 -= (m1 - m2) + (m3 + m4);
  m5 -= (m1 - m2) + (m3 - m4);
  m5 -= (m1 - m2) - (m3 + m4);
  m5 -= (m1 - m2) - (m3 - m4);
  UT_ASSERT( chkMatrixD22( m5, 2,4,6,8 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );

  m2 = 1*m1;
  m3 = m2*2;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 2,4,6,8 ) );

  m2 += 1*m1;
  m3 += m2*2;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 2,4,6,8 ) );
  UT_ASSERT( chkMatrixD22( m3, 6,12,18,24 ) );

  m2 -= 1*m1;
  m3 -= m2*2;
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 4,8,12,16 ) );

  m5 = 1*(m1 + m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 14,28,42,56 ) );
  m5 = 1*(m1 + m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 6,12,18,24 ) );
  m5 = 1*(m1 + m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -10,-20,-30,-40 ) );
  m5 = 1*(m1 + m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -2,-4,-6,-8 ) );
  m5 = 1*(m1 - m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 12,24,36,48 ) );
  m5 = 1*(m1 - m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  m5 = 1*(m1 - m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -12,-24,-36,-48 ) );
  m5 = 1*(m1 - m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -4,-8,-12,-16 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 4,8,12,16 ) );
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );

  m5 += 1*(m1 + m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 10,20,30,40 ) );
  m5 += 1*(m1 + m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 16,32,48,64 ) );
  m5 += 1*(m1 + m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 6,12,18,24 ) );
  m5 += 1*(m1 + m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  m5 += 1*(m1 - m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 16,32,48,64 ) );
  m5 += 1*(m1 - m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 20,40,60,80 ) );
  m5 += 1*(m1 - m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 8,16,24,32 ) );
  m5 += 1*(m1 - m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 4,8,12,16 ) );
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );

  m5 -= 1*(m1 + m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -10,-20,-30,-40 ) );
  m5 -= 1*(m1 + m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -16,-32,-48,-64 ) );
  m5 -= 1*(m1 + m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -6,-12,-18,-24 ) );
  m5 -= 1*(m1 + m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -4,-8,-12,-16 ) );
  m5 -= 1*(m1 - m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -16,-32,-48,-64 ) );
  m5 -= 1*(m1 - m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -20,-40,-60,-80 ) );
  m5 -= 1*(m1 - m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -8,-16,-24,-32 ) );
  m5 -= 1*(m1 - m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixD22( m5, -4,-8,-12,-16 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, 4,8,12,16 ) );
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );

  m5 = 1*(m1 + m2)*2;
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  m5 = 1*2*(m1 + m2);
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  m5 = (m1 + m2)*1*2;
  UT_ASSERT( chkMatrixD22( m5, 4,8,12,16 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );

  //Coverage of size function for s*(A+B)
  MatrixD<Real> m6 = 1*(m1 + m2)*2;
  UT_ASSERT( chkMatrixD22( m6, 4,8,12,16 ) );
  UT_ASSERT_EQUALS( m6.n(), 2 );
  UT_ASSERT_EQUALS( m6.m(), 2 );

  m2 = +m1;
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  m3 = -m2;
  UT_ASSERT( chkMatrixD22( m3, -1,-2,-3,-4 ) );
  m4 = +(m1 + m2);
  UT_ASSERT( chkMatrixD22( m4, 2,4,6,8 ) );
  m4 = +(m1 - m2);
  UT_ASSERT( chkMatrixD22( m4, 0,0,0,0 ) );
  m4 = -(m1 + m2);
  UT_ASSERT( chkMatrixD22( m4, -2,-4,-6,-8 ) );
  m4 = -(m1 - m2);
  UT_ASSERT( chkMatrixD22( m4, 0,0,0,0 ) );
  m4 = +(m1 + m2) + m3;
  UT_ASSERT( chkMatrixD22( m4, 1,2,3,4 ) );
  m4 = -(m1 + m2) + m3;
  UT_ASSERT( chkMatrixD22( m4, -3,-6,-9,-12 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m2, 1,2,3,4 ) );
  UT_ASSERT( chkMatrixD22( m3, -1,-2,-3,-4 ) );

  m4 = +1*m1;
  UT_ASSERT( chkMatrixD22( m4, 1,2,3,4 ) );
  m4 = -1*m1;
  UT_ASSERT( chkMatrixD22( m4, -1,-2,-3,-4 ) );
  m4 = +m1*1;
  UT_ASSERT( chkMatrixD22( m4, 1,2,3,4 ) );
  m4 = -m1*1;
  UT_ASSERT( chkMatrixD22( m4, -1,-2,-3,-4 ) );
  m4 = +(1*m1);
  UT_ASSERT( chkMatrixD22( m4, 1,2,3,4 ) );
  m4 = -(1*m1);
  UT_ASSERT( chkMatrixD22( m4, -1,-2,-3,-4 ) );
  m4 = +(m1*1);
  UT_ASSERT( chkMatrixD22( m4, 1,2,3,4 ) );
  m4 = -(m1*1);
  UT_ASSERT( chkMatrixD22( m4, -1,-2,-3,-4 ) );
  UT_ASSERT( chkMatrixD22( m1, 1,2,3,4 ) );
}
UT_TEST_CASE_END( MatrixD_ops2 )


//----------------------------------------------------------------------------//
// matrix-vector multiply
UT_TEST_CASE( MatrixD_Vector_Multiply2 )
{
  MatrixD<Real> m1(2, 2);
  MatrixD<Real> m2(2, 2);

  MatrixD<Real> col1(2, 1);
  MatrixD<Real> col2(2, 1);
  MatrixD<Real> col3(2, 1);

  MatrixD<Real> row1(1, 2);
  MatrixD<Real> row2(1, 2);
  MatrixD<Real> row3(1, 2);

  m1(0,0) = 3; m1(0,1) = 4;
  m1(1,0) = 5; m1(1,1) = 6;

  m2 = m1;

  col1(0,0) = 1;
  col1(1,0) = 2;

  col2 = col1;

  row1(0,0) = 1; row1(0,1) = 2;

  //Check column multiplication

  UT_ASSERT( chkMatrixD21( col1, 1,2 ) );
  UT_ASSERT( chkMatrixD21( col2, 1,2 ) );

  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );

  col2 = m1*col1;
  UT_ASSERT( chkMatrixD21( col2, 11,17 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );

  col2 += m1*col1;
  UT_ASSERT( chkMatrixD21( col2, 22,34 ) );

  col2 -= m1*col1;
  UT_ASSERT( chkMatrixD21( col2, 11,17 ) );


  col2 = (2*m1)*col1;
  UT_ASSERT( chkMatrixD21( col2, 2*11,2*17 ) );

  UT_ASSERT_EQUALS( ((2*m1)*col1).size(), 2 );

  col2 = +((2*m1)*col1);
  UT_ASSERT( chkMatrixD21( col2, 2*11,2*17 ) );

  col2 += (2*m1)*col1;
  UT_ASSERT( chkMatrixD21( col2, 4*11,4*17 ) );

  col2 = m1*(col1*2);
  UT_ASSERT( chkMatrixD21( col2, 2*11,2*17 ) );

  UT_ASSERT_EQUALS( (m1*(col1*2)).size(), 2 );

  col2 = +(m1*(col1*2));
  UT_ASSERT( chkMatrixD21( col2, 2*11,2*17 ) );

  col2 += m1*(col1*2);
  UT_ASSERT( chkMatrixD21( col2, 4*11,4*17 ) );



  col2 = (m1 + m2)*col1;
  UT_ASSERT( chkMatrixD21( col2, 22,34 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );

  col2 += (m1 + m2)*col1;
  UT_ASSERT( chkMatrixD21( col2, 2*22,2*34 ) );

  col2 -= (m1 + m2)*col1;
  UT_ASSERT( chkMatrixD21( col2, 22,34 ) );

  col2 = (m1 - m2)*col1;
  UT_ASSERT( chkMatrixD21( col2, 0,0 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );

  col2 = m1*col1;
  col2 += m1*col1 + m1*col1;
  UT_ASSERT( chkMatrixD21( col2, 3*11,3*17 ) );

  col2 += m1*col1 - m1*col1;
  UT_ASSERT( chkMatrixD21( col2, 3*11,3*17 ) );

  col2 = m1*col1 + col1;
  UT_ASSERT( chkMatrixD21( col2, 11+1,17+2 ) );


  m2(0,0) = 2; m2(0,1) = 3;
  m2(1,0) = 4; m2(1,1) = 5;

  col2(0,0) = 1;
  col2(1,0) = 2;

  col2 += (m1 - m2)*col1;
  UT_ASSERT( chkMatrixD21( col2, 4,5 ) );

  col2(0,0) = 7;
  col2(1,0) = 8;

  col2 -= (m1 - m2)*col1;
  UT_ASSERT( chkMatrixD21( col2, 4,5 ) );

  m2 = m1;
  col2 = col1;

  col3 = (m1 + m2)*(col1 + col2);
  UT_ASSERT( chkMatrixD21( col3, 44,68 ) );

  col3 += (m1 + m2)*(col1 + col2);
  UT_ASSERT( chkMatrixD21( col3, 2*44,2*68 ) );

  col3 -= (m1 + m2)*(col1 + col2);
  UT_ASSERT( chkMatrixD21( col3, 44,68 ) );

  col3 = (m1 + m2)*(col1 - col2);
  UT_ASSERT( chkMatrixD21( col3, 0,0 ) );

  col3 = (m1 - m2)*(col1 + col2);
  UT_ASSERT( chkMatrixD21( col3, 0,0 ) );

  col3 = (m1 - m2)*(col1 - col2);
  UT_ASSERT( chkMatrixD21( col3, 0,0 ) );

  MatrixD<Real> col4 = (m1 - m2)*(col1 - col2);
  UT_ASSERT( chkMatrixD21( col4, 0,0 ) );


  //Just to get complete coverage with the unary operator+
  col2 = +(m1*col1);
  UT_ASSERT( chkMatrixD21( col2, 11,17 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );

  col2 = +((m1 + m2)*col1);
  UT_ASSERT( chkMatrixD21( col2, 22,34 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );



  //Check row multiplication
  row2 = row1*m1;
  UT_ASSERT( chkMatrixD12( row2, 13,16 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );

  row2 = (2*row1)*m1;
  UT_ASSERT( chkMatrixD12( row2, 2*13,2*16 ) );

  row2 = row1*(m1*2);
  UT_ASSERT( chkMatrixD12( row2, 2*13,2*16 ) );

  row2 = row1*(m1 + m2);
  UT_ASSERT( chkMatrixD12( row2, 26,32 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );

  row2 += row1*(m1 + m2);
  UT_ASSERT( chkMatrixD12( row2, 2*26,2*32 ) );

  row2 -= row1*(m1 + m2);
  UT_ASSERT( chkMatrixD12( row2, 26,32 ) );

  row2 = row1*(m1 - m2);
  UT_ASSERT( chkMatrixD12( row2, 0,0 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );

  //Just to get complete coverage with the unary operator+
  row2 = +(row1*(m1 + m2));
  UT_ASSERT( chkMatrixD12( row2, 26,32 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );

  //Coverage of size function in A*(B + C)
  MatrixD<Real> row5 = row1*(m1 + m2);
  UT_ASSERT( chkMatrixD12( row5, 26,32 ) );
  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );
  UT_ASSERT_EQUALS( row5.m(), 1 );
  UT_ASSERT_EQUALS( row5.n(), 2 );



  row2 = row1;

  row3 = (row1 + row2)*(m1 + m2);
  UT_ASSERT( chkMatrixD12( row3, 52,64 ) );

  row3 = (row1 - row2)*(m1 + m2);
  UT_ASSERT( chkMatrixD12( row3, 0,0 ) );

  row3 = (row1 + row2)*(m1 - m2);
  UT_ASSERT( chkMatrixD12( row3, 0,0 ) );

  row3 = (row1 - row2)*(m1 - m2);
  UT_ASSERT( chkMatrixD12( row3, 0,0 ) );
  UT_ASSERT( chkMatrixD12( row1, 1,2 ) );
  UT_ASSERT( chkMatrixD12( row2, 1,2 ) );

  UT_ASSERT( chkMatrixD22( m1, 3,4,5,6 ) );
  UT_ASSERT( chkMatrixD22( m2, 3,4,5,6 ) );

  //Just to get complete coverage with the unary operator+
  row3 = +( (row1 + row2)*(m1 + m2) );
  UT_ASSERT( chkMatrixD12( row3, 52,64 ) );


  MatrixD<Real> val(1, 1);

  val = row1*col1;
  UT_ASSERT_EQUALS( val(0,0), 5 );

  val = row1*m1*col1;
  UT_ASSERT_EQUALS( val(0,0), 45 );


  //It is ok for what is on the left to be on the right as long as it is the first thing and not part of a multiplication
  col2 = col1;
  col2 = col2 + m1*col1;
  UT_ASSERT( chkMatrixD21( col2, 11+1,17+2 ) );

  //This is NOT ok because of lazy expressions the variable left of '=' cannot be part of a multiplication on the right of '='
  //UT_ASSERT_THROW( col2 = col1 + m1*col2;, AssertionException );
  //UT_ASSERT_THROW( col2 = m1*col1 + col2;, AssertionException );

  //Make sure that an exception is thrown when dimensions do not match
  //or when the same veriable is used on the left and right as part of a multiplication

  //UT_ASSERT_THROW( row2  = m1*col1, AssertionException );
  //UT_ASSERT_THROW( row2 += m1*col1, AssertionException );
  //UT_ASSERT_THROW( row2 -= m1*col1, AssertionException );

  //UT_ASSERT_THROW( col2  = row1*m1, AssertionException );
  //UT_ASSERT_THROW( col2 += row1*m1, AssertionException );
  //UT_ASSERT_THROW( col2 -= row1*m1, AssertionException );

  //UT_ASSERT_THROW( row2  = col1*m1, AssertionException );
  //UT_ASSERT_THROW( row2 += col1*m1, AssertionException );
  //UT_ASSERT_THROW( row2 -= col1*m1, AssertionException );

  //UT_ASSERT_THROW( col2  = m1*row1, AssertionException );
  //UT_ASSERT_THROW( col2 += m1*row1, AssertionException );
  //UT_ASSERT_THROW( col2 -= m1*row1, AssertionException );

  //UT_ASSERT_THROW( m2  = m2*m1, AssertionException );
  //UT_ASSERT_THROW( m2 += m2*m1, AssertionException );
  //UT_ASSERT_THROW( m2 -= m2*m1, AssertionException );

  //UT_ASSERT_THROW( m2  = m1*m2, AssertionException );
  //UT_ASSERT_THROW( m2 += m1*m2, AssertionException );
  //UT_ASSERT_THROW( m2 -= m1*m2, AssertionException );

  //UT_ASSERT_THROW( col1  = m1*col1, AssertionException );
  //UT_ASSERT_THROW( col1 += m1*col1, AssertionException );
  //UT_ASSERT_THROW( col1 -= m1*col1, AssertionException );

  //UT_ASSERT_THROW( row1  = row1*m1, AssertionException );
  //UT_ASSERT_THROW( row1 += row1*m1, AssertionException );
  //UT_ASSERT_THROW( row1 -= row1*m1, AssertionException );

}
UT_TEST_CASE_END( MatrixD_Vector_Multiply2 )


//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixD_sub_matrix )
{
  {
  MatrixD<Real> M1(2, 2);
  M1 = 0;
  M1(0,0) = 1;
  M1(1,0) = 3;
  M1(1,1) = 2;

  MatrixDView<Real> row0 = M1.row(0);
  MatrixDView<Real> row1 = M1.row(1);

  UT_ASSERT_EQUALS( 2, row0.size() );
  UT_ASSERT_EQUALS( 1, row0.m() );
  UT_ASSERT_EQUALS( 2, row0.n() );
  UT_ASSERT_EQUALS( 2, row0.stride() );

  UT_ASSERT_EQUALS( 2, row1.size() );
  UT_ASSERT_EQUALS( 1, row1.m() );
  UT_ASSERT_EQUALS( 2, row1.n() );
  UT_ASSERT_EQUALS( 2, row1.stride() );

  UT_ASSERT_EQUALS( 1, row0(0,0) );
  UT_ASSERT_EQUALS( 0, row0(0,1) );
  UT_ASSERT_EQUALS( 3, row1(0,0) );
  UT_ASSERT_EQUALS( 2, row1(0,1) );

  MatrixDView<Real> col0 = M1.col(0);
  MatrixDView<Real> col1 = M1.col(1);

  UT_ASSERT_EQUALS( 2, col0.size() );
  UT_ASSERT_EQUALS( 2, col0.m() );
  UT_ASSERT_EQUALS( 1, col0.n() );
  UT_ASSERT_EQUALS( 2, col0.stride() );

  UT_ASSERT_EQUALS( 2, col1.size() );
  UT_ASSERT_EQUALS( 2, col1.m() );
  UT_ASSERT_EQUALS( 1, col1.n() );
  UT_ASSERT_EQUALS( 2, col1.stride() );

  UT_ASSERT_EQUALS( 1, col0(0,0) );
  UT_ASSERT_EQUALS( 0, col1(0,0) );
  UT_ASSERT_EQUALS( 3, col0(1,0) );
  UT_ASSERT_EQUALS( 2, col1(1,0) );

  MatrixDView<Real> sub = M1.sub(1,1,1,1);

  UT_ASSERT_EQUALS( 1, sub.size() );
  UT_ASSERT_EQUALS( 1, sub.m() );
  UT_ASSERT_EQUALS( 1, sub.n() );
  UT_ASSERT_EQUALS( 2, col0.stride() );

  UT_ASSERT_EQUALS( 2, sub(0,0) );
  }

  {
  MatrixD<Real> M2(4,2);
  M2 = 0;
  M2(0,0) = 2, M2(0,1) = 9;
  M2(1,0) = 3, M2(1,1) = 8;
  M2(2,0) = 4, M2(2,1) = 7;
  M2(3,0) = 5, M2(3,1) = 6;

  MatrixDView<Real> row0 = M2.row(0);
  MatrixDView<Real> row1 = M2.row(1);
  MatrixDView<Real> row2 = M2.row(2);
  MatrixDView<Real> row3 = M2.row(3);

  UT_ASSERT_EQUALS( 2, row0.size() );
  UT_ASSERT_EQUALS( 1, row0.m() );
  UT_ASSERT_EQUALS( 2, row0.n() );
  UT_ASSERT_EQUALS( 2, row0.stride() );

  UT_ASSERT_EQUALS( 2, row1.size() );
  UT_ASSERT_EQUALS( 1, row1.m() );
  UT_ASSERT_EQUALS( 2, row1.n() );
  UT_ASSERT_EQUALS( 2, row1.stride() );

  UT_ASSERT_EQUALS( 2, row2.size() );
  UT_ASSERT_EQUALS( 1, row2.m() );
  UT_ASSERT_EQUALS( 2, row2.n() );
  UT_ASSERT_EQUALS( 2, row2.stride() );

  UT_ASSERT_EQUALS( 2, row3.size() );
  UT_ASSERT_EQUALS( 1, row3.m() );
  UT_ASSERT_EQUALS( 2, row3.n() );
  UT_ASSERT_EQUALS( 2, row3.stride() );

  UT_ASSERT_EQUALS( 2, row0(0,0) );
  UT_ASSERT_EQUALS( 9, row0(0,1) );
  UT_ASSERT_EQUALS( 3, row1(0,0) );
  UT_ASSERT_EQUALS( 8, row1(0,1) );
  UT_ASSERT_EQUALS( 4, row2(0,0) );
  UT_ASSERT_EQUALS( 7, row2(0,1) );
  UT_ASSERT_EQUALS( 5, row3(0,0) );
  UT_ASSERT_EQUALS( 6, row3(0,1) );

  MatrixDView<Real> col0 = M2.col(0);
  MatrixDView<Real> col1 = M2.col(1);

  UT_ASSERT_EQUALS( 4, col0.size() );
  UT_ASSERT_EQUALS( 4, col0.m() );
  UT_ASSERT_EQUALS( 1, col0.n() );
  UT_ASSERT_EQUALS( 4, col0.size() );

  UT_ASSERT_EQUALS( 4, col1.size() );
  UT_ASSERT_EQUALS( 4, col1.m() );
  UT_ASSERT_EQUALS( 1, col1.n() );
  UT_ASSERT_EQUALS( 4, col1.size() );

  UT_ASSERT_EQUALS( 2, col0(0,0) );
  UT_ASSERT_EQUALS( 3, col0(1,0) );
  UT_ASSERT_EQUALS( 4, col0(2,0) );
  UT_ASSERT_EQUALS( 5, col0(3,0) );
  UT_ASSERT_EQUALS( 9, col1(0,0) );
  UT_ASSERT_EQUALS( 8, col1(1,0) );
  UT_ASSERT_EQUALS( 7, col1(2,0) );
  UT_ASSERT_EQUALS( 6, col1(3,0) );
  }

  {
  MatrixD<Real> M3(2,4);
  M3 = 0;
  M3(0,0) = 2, M3(0,1) = 3, M3(0,2) = 4, M3(0,3) = 5;
  M3(1,0) = 9, M3(1,1) = 8, M3(1,2) = 7, M3(1,3) = 6;

  MatrixDView<Real> row0 = M3.row(0);
  MatrixDView<Real> row1 = M3.row(1);

  UT_ASSERT_EQUALS( 4, row0.size() );
  UT_ASSERT_EQUALS( 1, row0.m() );
  UT_ASSERT_EQUALS( 4, row0.n() );
  UT_ASSERT_EQUALS( 4, row0.stride() );

  UT_ASSERT_EQUALS( 4, row1.size() );
  UT_ASSERT_EQUALS( 1, row1.m() );
  UT_ASSERT_EQUALS( 4, row1.n() );
  UT_ASSERT_EQUALS( 4, row1.stride() );

  UT_ASSERT_EQUALS( 2, row0(0,0) );
  UT_ASSERT_EQUALS( 3, row0(0,1) );
  UT_ASSERT_EQUALS( 4, row0(0,2) );
  UT_ASSERT_EQUALS( 5, row0(0,3) );

  UT_ASSERT_EQUALS( 9, row1(0,0) );
  UT_ASSERT_EQUALS( 8, row1(0,1) );
  UT_ASSERT_EQUALS( 7, row1(0,2) );
  UT_ASSERT_EQUALS( 6, row1(0,3) );

  MatrixDView<Real> col0 = M3.col(0);
  MatrixDView<Real> col1 = M3.col(1);
  MatrixDView<Real> col2 = M3.col(2);
  MatrixDView<Real> col3 = M3.col(3);

  UT_ASSERT_EQUALS( 2, col0.size() );
  UT_ASSERT_EQUALS( 2, col0.m() );
  UT_ASSERT_EQUALS( 1, col0.n() );
  UT_ASSERT_EQUALS( 2, col0.size() );

  UT_ASSERT_EQUALS( 2, col1.size() );
  UT_ASSERT_EQUALS( 2, col1.m() );
  UT_ASSERT_EQUALS( 1, col1.n() );
  UT_ASSERT_EQUALS( 2, col1.size() );

  UT_ASSERT_EQUALS( 2, col2.size() );
  UT_ASSERT_EQUALS( 2, col2.m() );
  UT_ASSERT_EQUALS( 1, col2.n() );
  UT_ASSERT_EQUALS( 2, col2.size() );

  UT_ASSERT_EQUALS( 2, col3.size() );
  UT_ASSERT_EQUALS( 2, col3.m() );
  UT_ASSERT_EQUALS( 1, col3.n() );
  UT_ASSERT_EQUALS( 2, col3.size() );

  UT_ASSERT_EQUALS( 2, col0(0,0) );
  UT_ASSERT_EQUALS( 9, col0(1,0) );
  UT_ASSERT_EQUALS( 3, col1(0,0) );
  UT_ASSERT_EQUALS( 8, col1(1,0) );
  UT_ASSERT_EQUALS( 4, col2(0,0) );
  UT_ASSERT_EQUALS( 7, col2(1,0) );
  UT_ASSERT_EQUALS( 5, col3(0,0) );
  UT_ASSERT_EQUALS( 6, col3(1,0) );
  }

  {
  MatrixD<Real> M3(3,4);
  M3 = 0;
  M3(0,0) =  2, M3(0,1) =  3, M3(0,2) =  4, M3(0,3) =  5;
  M3(1,0) =  9, M3(1,1) =  8, M3(1,2) =  7, M3(1,3) =  6;
  M3(2,0) = 10, M3(2,1) = 11, M3(2,2) = 12, M3(2,3) = 13;

  MatrixD<Real> sub(2,3);

  subMatrixValue(M3, {0,2}, {0,1,3}, 2., sub);

  UT_ASSERT_EQUALS( 2*2, sub(0,0) );
  UT_ASSERT_EQUALS( 2*3, sub(0,1) );
  UT_ASSERT_EQUALS( 2*5, sub(0,2) );

  UT_ASSERT_EQUALS( 2*10, sub(1,0) );
  UT_ASSERT_EQUALS( 2*11, sub(1,1) );
  UT_ASSERT_EQUALS( 2*13, sub(1,2) );

  subMatrixPlus(M3, {0,2}, {0,1,3}, 3., sub);

  UT_ASSERT_EQUALS( 3*2+2*2, sub(0,0) );
  UT_ASSERT_EQUALS( 3*3+2*3, sub(0,1) );
  UT_ASSERT_EQUALS( 3*5+2*5, sub(0,2) );

  UT_ASSERT_EQUALS( 3*10+2*10, sub(1,0) );
  UT_ASSERT_EQUALS( 3*11+2*11, sub(1,1) );
  UT_ASSERT_EQUALS( 3*13+2*13, sub(1,2) );
  }
}
UT_TEST_CASE_END( MatrixD_sub_matrix )


//############################################################################//
UT_TEST_SUITE_END(DenseLinAlg_MatrixD)
