// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// MatrixS_btest
// testing of MatrixS<M,N,T> class

#include "unit_tester.hpp"

#include "numpack/DenseLinAlg/DynamicSize/MatrixD.h"
#include "numpack/DenseLinAlg/DynamicSize/MatrixD_Sub.h"
#include "numpack/DenseLinAlg/DynamicSize/MatrixSymD.h"
#include "numpack/DenseLinAlg/DynamicSize/VectorD.h"

#include "chkMatrixD_btest.h"

#include <iostream>

namespace numpack
{
namespace DLA
{

//template class MatrixSymD<int>;
//template class MatrixSymD<int>;

/*
typedef MatrixS<1,1,int> Mat;
typedef OpAddS< Mat, Mat, true > MatAddT;

template class OpMulS< Mat, Mat >;
template class OpMulS< MatAddT, Mat >;
template class OpMulS< Mat, MatAddT >;
template class OpMulS< MatAddT, MatAddT >;

template class OpAddS< Mat, Mat, true >;
template class OpSubS< Mat, Mat, true >;
template class OpMulSScalar< Mat, Real, true >;
template class OpMulSFactor< Mat, OpMulSScalar< Mat, Real, true > >;
template class OpMulSFactor< OpMulSScalar< Mat, Real, true >, Mat >;

typedef OpAddS< Mat, Mat, false > MatAddF;

template class OpAddS< Mat, Mat, false >;
template class OpSubS< Mat, Mat, false >;
template class OpMulSScalar< Mat, Real, false >;
template class OpMulSFactor< Mat, OpMulSScalar< Mat, Real, false > >;
template class OpMulSFactor< OpMulSScalar< Mat, Real, false >, Mat >;
*/
}
}

using namespace numpack;
using namespace numpack::DLA;

//############################################################################//
UT_TEST_SUITE( MatrixSymD_tests )

typedef int Int;

typedef VectorD<Int> VectorD;
typedef MatrixD<Int> MatrixD;
typedef MatrixSymD<Int> MatrixSymD;


UT_TEST_CASE( test1 )
{
  const int n = 2;

  DLA::MatrixSymD<Real> m1(n);
  DLA::MatrixSymD<Real> m2(n);

  m1 = 1;
  m2 = 2;

  DLA::MatrixSymD<Real> m3(n);
  m3 =  m1+m2;
  m3.dump();

  DLA::MatrixD<Real> m4(2,2);

  //m4 = m1*m2;
  //m4.dump();


}
UT_TEST_CASE_END( test1 )

#if 0
//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( size1 )
{
  BOOST_CHECK( MatrixSymS1::M == 1 );
  BOOST_CHECK( MatrixSymS1::N == 1 );
}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( matrix_ops1 )
{
  const Int data = 3;
  MatrixSymS1 m1(&data, 1);
  MatrixSymS1 m2(m1);
  MatrixSymS1 m3(1), m4{2}, m5;

  // size
  BOOST_CHECK( m1.M == 1 );
  BOOST_CHECK( m1.N == 1 );
  BOOST_CHECK( m2.M == 1 );
  BOOST_CHECK( m2.N == 1 );

  // ctors
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  1, m3(0,0) );
  BOOST_CHECK_EQUAL(  2, m4(0,0) );

  // assignment
  m3 = m1;
  m4 = data;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m3(0,0) );
  BOOST_CHECK_EQUAL(  3, m4(0,0) );

  m1 = m2 = m3 = data;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  3, m3(0,0) );

  m4 = data;
  m1 = m2 = m3 = m4;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  3, m3(0,0) );
  BOOST_CHECK_EQUAL(  3, m4(0,0) );

  // unary
  m2 = +m1;
  m3 = -m1;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL( -3, m3(0,0) );

  // binary accumulation
  m3 = m1;
  m4 = m1;
  m3 *= data;
  m4 /= data;
  BOOST_CHECK_EQUAL(  9, m3(0,0) );
  BOOST_CHECK_EQUAL(  1, m4(0,0) );

  m1 = data;
  m2 = m1;
  m3 = m1;
  m2 += m1;
  m3 -= m1;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  6, m2(0,0) );
  BOOST_CHECK_EQUAL(  0, m3(0,0) );

  // binary operators
  m1 = data;
//  m2 = m1 + data;
//  m3 = m1 - data;
  m4 = m1 * data;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
//  BOOST_CHECK_EQUAL(  6, m2(0,0) );
//  BOOST_CHECK_EQUAL(  0, m3(0,0) );
  BOOST_CHECK_EQUAL(  9, m4(0,0) );

//  m2 = data + m1;
//  m3 = data - m1;
  m4 = data * m1;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
//  BOOST_CHECK_EQUAL(  6, m2(0,0) );
//  BOOST_CHECK_EQUAL(  0, m3(0,0) );
  BOOST_CHECK_EQUAL(  9, m4(0,0) );

  m1 = m2 = data;
  m3 = m1 + m2;
  m4 = m1 - m2;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  6, m3(0,0) );
  BOOST_CHECK_EQUAL(  0, m4(0,0) );

  // arithmetic combinations

  m1 = m2 = data;
  m3 = m1 + m2;
  m4 = m1 + m2 + m3;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  6, m3(0,0) );
  BOOST_CHECK_EQUAL( 12, m4(0,0) );

  m2 += m1;
  m3 += m1 + m2;
  m4 += m1 + m2 + m3;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  6, m2(0,0) );
  BOOST_CHECK_EQUAL( 15, m3(0,0) );
  BOOST_CHECK_EQUAL( 36, m4(0,0) );

  m3 = m1 - m2;
  m4 = m1 - m2 - m3;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  6, m2(0,0) );
  BOOST_CHECK_EQUAL( -3, m3(0,0) );
  BOOST_CHECK_EQUAL(  0, m4(0,0) );

  m2 -= m1;
  m3 -= m1 - m2;
  m4 -= m1 - m2 - m3;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL( -3, m3(0,0) );
  BOOST_CHECK_EQUAL( -3, m4(0,0) );

  m3 = m1 - m2;
  m4 = m1 + m2 - m3;
  m5 = m1 - m2 + m3;
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  0, m3(0,0) );
  BOOST_CHECK_EQUAL(  6, m4(0,0) );
  BOOST_CHECK_EQUAL(  0, m5(0,0) );

  m5 = (m1 + m2) + (m3 + m4);
  BOOST_CHECK_EQUAL( 12, m5(0,0) );
  m5 = (m1 + m2) + (m3 - m4);
  BOOST_CHECK_EQUAL(  0, m5(0,0) );
  m5 = (m1 + m2) - (m3 + m4);
  BOOST_CHECK_EQUAL(  0, m5(0,0) );
  m5 = (m1 + m2) - (m3 - m4);
  BOOST_CHECK_EQUAL( 12, m5(0,0) );
  m5 = (m1 - m2) + (m3 + m4);
  BOOST_CHECK_EQUAL(  6, m5(0,0) );
  m5 = (m1 - m2) + (m3 - m4);
  BOOST_CHECK_EQUAL( -6, m5(0,0) );
  m5 = (m1 - m2) - (m3 + m4);
  BOOST_CHECK_EQUAL( -6, m5(0,0) );
  m5 = (m1 - m2) - (m3 - m4);
  BOOST_CHECK_EQUAL(  6, m5(0,0) );
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  0, m3(0,0) );
  BOOST_CHECK_EQUAL(  6, m4(0,0) );

  m5 += (m1 + m2) + (m3 + m4);
  m5 += (m1 + m2) + (m3 - m4);
  m5 += (m1 + m2) - (m3 + m4);
  m5 += (m1 + m2) - (m3 - m4);
  m5 += (m1 - m2) + (m3 + m4);
  m5 += (m1 - m2) + (m3 - m4);
  m5 += (m1 - m2) - (m3 + m4);
  m5 += (m1 - m2) - (m3 - m4);
  BOOST_CHECK_EQUAL( 30, m5(0,0) );
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  0, m3(0,0) );
  BOOST_CHECK_EQUAL(  6, m4(0,0) );

  m5 -= (m1 + m2) + (m3 + m4);
  m5 -= (m1 + m2) + (m3 - m4);
  m5 -= (m1 + m2) - (m3 + m4);
  m5 -= (m1 + m2) - (m3 - m4);
  m5 -= (m1 - m2) + (m3 + m4);
  m5 -= (m1 - m2) + (m3 - m4);
  m5 -= (m1 - m2) - (m3 + m4);
  m5 -= (m1 - m2) - (m3 - m4);
  BOOST_CHECK_EQUAL(  6, m5(0,0) );
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL(  0, m3(0,0) );
  BOOST_CHECK_EQUAL(  6, m4(0,0) );

  m1 = data;

  m2 = 4*m1;
  m3 = m2*7;
  BOOST_CHECK_EQUAL(   3, m1(0,0) );
  BOOST_CHECK_EQUAL(  12, m2(0,0) );
  BOOST_CHECK_EQUAL(  84, m3(0,0) );

  m2 += 4*m1;
  m3 += m2*7;
  BOOST_CHECK_EQUAL(   3, m1(0,0) );
  BOOST_CHECK_EQUAL(  24, m2(0,0) );
  BOOST_CHECK_EQUAL( 252, m3(0,0) );

  m2 -= 4*m1;
  m3 -= m2*7;
  BOOST_CHECK_EQUAL(   3, m1(0,0) );
  BOOST_CHECK_EQUAL(  12, m2(0,0) );
  BOOST_CHECK_EQUAL( 168, m3(0,0) );

  m5 = 2*(m1 + m2) + (m3 + m4)*3;
  BOOST_CHECK_EQUAL(  552, m5(0,0) );
  m5 = 2*(m1 + m2) + (m3 - m4)*3;
  BOOST_CHECK_EQUAL(  516, m5(0,0) );
  m5 = 2*(m1 + m2) - (m3 + m4)*3;
  BOOST_CHECK_EQUAL( -492, m5(0,0) );
  m5 = 2*(m1 + m2) - (m3 - m4)*3;
  BOOST_CHECK_EQUAL( -456, m5(0,0) );
  m5 = 2*(m1 - m2) + (m3 + m4)*3;
  BOOST_CHECK_EQUAL(  504, m5(0,0) );
  m5 = 2*(m1 - m2) + (m3 - m4)*3;
  BOOST_CHECK_EQUAL(  468, m5(0,0) );
  m5 = 2*(m1 - m2) - (m3 + m4)*3;
  BOOST_CHECK_EQUAL( -540, m5(0,0) );
  m5 = 2*(m1 - m2) - (m3 - m4)*3;
  BOOST_CHECK_EQUAL( -504, m5(0,0) );
  BOOST_CHECK_EQUAL(    3, m1(0,0) );
  BOOST_CHECK_EQUAL(   12, m2(0,0) );
  BOOST_CHECK_EQUAL(  168, m3(0,0) );
  BOOST_CHECK_EQUAL(    6, m4(0,0) );

  m5 += 2*(m1 + m2) + (m3 + m4)*3;
  BOOST_CHECK_EQUAL(   48, m5(0,0) );
  m5 += 2*(m1 + m2) + (m3 - m4)*3;
  BOOST_CHECK_EQUAL(  564, m5(0,0) );
  m5 += 2*(m1 + m2) - (m3 + m4)*3;
  BOOST_CHECK_EQUAL(   72, m5(0,0) );
  m5 += 2*(m1 + m2) - (m3 - m4)*3;
  BOOST_CHECK_EQUAL( -384, m5(0,0) );
  m5 += 2*(m1 - m2) + (m3 + m4)*3;
  BOOST_CHECK_EQUAL(  120, m5(0,0) );
  m5 += 2*(m1 - m2) + (m3 - m4)*3;
  BOOST_CHECK_EQUAL(  588, m5(0,0) );
  m5 += 2*(m1 - m2) - (m3 + m4)*3;
  BOOST_CHECK_EQUAL(   48, m5(0,0) );
  m5 += 2*(m1 - m2) - (m3 - m4)*3;
  BOOST_CHECK_EQUAL( -456, m5(0,0) );
  BOOST_CHECK_EQUAL(    3, m1(0,0) );
  BOOST_CHECK_EQUAL(   12, m2(0,0) );
  BOOST_CHECK_EQUAL(  168, m3(0,0) );
  BOOST_CHECK_EQUAL(    6, m4(0,0) );

  m5 -= 2*(m1 + m2) + (m3 + m4)*3;
  BOOST_CHECK_EQUAL( -1008, m5(0,0) );
  m5 -= 2*(m1 + m2) + (m3 - m4)*3;
  BOOST_CHECK_EQUAL( -1524, m5(0,0) );
  m5 -= 2*(m1 + m2) - (m3 + m4)*3;
  BOOST_CHECK_EQUAL( -1032, m5(0,0) );
  m5 -= 2*(m1 + m2) - (m3 - m4)*3;
  BOOST_CHECK_EQUAL(  -576, m5(0,0) );
  m5 -= 2*(m1 - m2) + (m3 + m4)*3;
  BOOST_CHECK_EQUAL( -1080, m5(0,0) );
  m5 -= 2*(m1 - m2) + (m3 - m4)*3;
  BOOST_CHECK_EQUAL( -1548, m5(0,0) );
  m5 -= 2*(m1 - m2) - (m3 + m4)*3;
  BOOST_CHECK_EQUAL( -1008, m5(0,0) );
  m5 -= 2*(m1 - m2) - (m3 - m4)*3;
  BOOST_CHECK_EQUAL(  -504, m5(0,0) );
  BOOST_CHECK_EQUAL(     3, m1(0,0) );
  BOOST_CHECK_EQUAL(    12, m2(0,0) );
  BOOST_CHECK_EQUAL(   168, m3(0,0) );
  BOOST_CHECK_EQUAL(     6, m4(0,0) );

  m5 = 2*(m1 + m2)*3;
  BOOST_CHECK_EQUAL( 90, m5(0,0) );
  m5 = 2*3*(m1 + m2);
  BOOST_CHECK_EQUAL( 90, m5(0,0) );
  m5 = (m1 + m2)*2*3;
  BOOST_CHECK_EQUAL( 90, m5(0,0) );
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL( 12, m2(0,0) );

  m2 = +m1;
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  m3 = -m2;
  BOOST_CHECK_EQUAL( -3, m3(0,0) );
  m4 = +(m1 + m2);
  BOOST_CHECK_EQUAL(  6, m4(0,0) );
  m4 = +(m1 - m2);
  BOOST_CHECK_EQUAL(  0, m4(0,0) );
  m4 = -(m1 + m2);
  BOOST_CHECK_EQUAL( -6, m4(0,0) );
  m4 = -(m1 - m2);
  BOOST_CHECK_EQUAL(  0, m4(0,0) );
  m4 = +(m1 + m2) + m3;
  BOOST_CHECK_EQUAL(  3, m4(0,0) );
  m4 = -(m1 + m2) + m3;
  BOOST_CHECK_EQUAL( -9, m4(0,0) );
  BOOST_CHECK_EQUAL(  3, m1(0,0) );
  BOOST_CHECK_EQUAL(  3, m2(0,0) );
  BOOST_CHECK_EQUAL( -3, m3(0,0) );

  m4 = +5*m1;
  BOOST_CHECK_EQUAL(  15, m4(0,0) );
  m4 = -5*m1;
  BOOST_CHECK_EQUAL( -15, m4(0,0) );
  m4 = +m1*5;
  BOOST_CHECK_EQUAL(  15, m4(0,0) );
  m4 = -m1*5;
  BOOST_CHECK_EQUAL( -15, m4(0,0) );
  m4 = +(5*m1);
  BOOST_CHECK_EQUAL(  15, m4(0,0) );
  m4 = -(5*m1);
  BOOST_CHECK_EQUAL( -15, m4(0,0) );
  m4 = +(m1*5);
  BOOST_CHECK_EQUAL(  15, m4(0,0) );
  m4 = -(m1*5);
  BOOST_CHECK_EQUAL( -15, m4(0,0) );
  BOOST_CHECK_EQUAL(   3, m1(0,0) );

  m4 = {6};
  BOOST_CHECK_EQUAL( 6, m4(0,0) );

  Int i0 = m1*m1;
  BOOST_CHECK_EQUAL( 9, i0 );

  Int i1 = (m1+m1)*m1;
  BOOST_CHECK_EQUAL( 18, i1 );

  Int i2 = m1*(m1+m1);
  BOOST_CHECK_EQUAL( 18, i2 );

  Int i3 = (m1+m1)*(m1+m1);
  BOOST_CHECK_EQUAL( 36, i3 );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( size2 )
{
  BOOST_CHECK( MatrixSymS2::M == 2 );
  BOOST_CHECK( MatrixSymS2::N == 2 );
}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( identity2 )
{
  MatrixSymS2 m1;
  m1 = Identity();

  BOOST_CHECK_EQUAL(  1, m1(0,0) );
  BOOST_CHECK_EQUAL(  0, m1(0,1) );
  BOOST_CHECK_EQUAL(  0, m1(1,0) );
  BOOST_CHECK_EQUAL(  1, m1(1,1) );

  MatrixSymS2 m2;
  MatrixSymS2 I = Identity();

  m2 = 2*I;

  BOOST_CHECK_EQUAL(  2, m2(0,0) );
  BOOST_CHECK_EQUAL(  0, m2(0,1) );
  BOOST_CHECK_EQUAL(  0, m2(1,0) );
  BOOST_CHECK_EQUAL(  2, m2(1,1) );

  MatrixSymS2 m3 = Identity();

  BOOST_CHECK_EQUAL(  1, m3(0,0) );
  BOOST_CHECK_EQUAL(  0, m3(0,1) );
  BOOST_CHECK_EQUAL(  0, m3(1,0) );
  BOOST_CHECK_EQUAL(  1, m3(1,1) );

}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( MatrixSymS_3x3_test )
{
  Int a = 1; Int b = 4; Int c = 7;
  Int d = 4; Int e = 5; Int f = 8;
  Int g = 7; Int h = 8; Int i = 9;

  Int Adata[] = {a,
                 d, e,
                 g, h, i};

  MatrixSymS<3,Int> A(Adata, 6);

  BOOST_CHECK_EQUAL(A(0,0), a);
  BOOST_CHECK_EQUAL(A(0,1), b);
  BOOST_CHECK_EQUAL(A(0,2), c);

  BOOST_CHECK_EQUAL(A(1,0), d);
  BOOST_CHECK_EQUAL(A(1,1), e);
  BOOST_CHECK_EQUAL(A(1,2), f);

  BOOST_CHECK_EQUAL(A(2,0), g);
  BOOST_CHECK_EQUAL(A(2,1), h);
  BOOST_CHECK_EQUAL(A(2,2), i);
}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( ATransposeA_ctor )
{
  const Int data[4] = {1,2,3,4};
  MatrixS2 m1(data, 4), m2(m1);

  MatrixSymS2 sm1a;
  sm1a = Transpose(m1)*m1;

  BOOST_CHECK_EQUAL(  10, sm1a(0,0) );
  BOOST_CHECK_EQUAL(  14, sm1a(0,1) );
  BOOST_CHECK_EQUAL(  14, sm1a(1,0) );
  BOOST_CHECK_EQUAL(  20, sm1a(1,1) );

  sm1a += Transpose(m1)*m1;

  BOOST_CHECK_EQUAL( 2*10, sm1a(0,0) );
  BOOST_CHECK_EQUAL( 2*14, sm1a(0,1) );
  BOOST_CHECK_EQUAL( 2*14, sm1a(1,0) );
  BOOST_CHECK_EQUAL( 2*20, sm1a(1,1) );

  sm1a -= Transpose(m1)*m1;

  BOOST_CHECK_EQUAL(  10, sm1a(0,0) );
  BOOST_CHECK_EQUAL(  14, sm1a(0,1) );
  BOOST_CHECK_EQUAL(  14, sm1a(1,0) );
  BOOST_CHECK_EQUAL(  20, sm1a(1,1) );

  MatrixSymS2 sm1b;
  sm1b = m1*Transpose(m1);

  BOOST_CHECK_EQUAL(   5, sm1b(0,0) );
  BOOST_CHECK_EQUAL(  11, sm1b(0,1) );
  BOOST_CHECK_EQUAL(  11, sm1b(1,0) );
  BOOST_CHECK_EQUAL(  25, sm1b(1,1) );

  sm1b += m1*Transpose(m1);

  BOOST_CHECK_EQUAL(  2*5, sm1b(0,0) );
  BOOST_CHECK_EQUAL( 2*11, sm1b(0,1) );
  BOOST_CHECK_EQUAL( 2*11, sm1b(1,0) );
  BOOST_CHECK_EQUAL( 2*25, sm1b(1,1) );

  sm1b -= m1*Transpose(m1);

  BOOST_CHECK_EQUAL(   5, sm1b(0,0) );
  BOOST_CHECK_EQUAL(  11, sm1b(0,1) );
  BOOST_CHECK_EQUAL(  11, sm1b(1,0) );
  BOOST_CHECK_EQUAL(  25, sm1b(1,1) );


  MatrixSymS2 sm2a = Transpose(m1)*m1;

  BOOST_CHECK_EQUAL(  10, sm2a(0,0) );
  BOOST_CHECK_EQUAL(  14, sm2a(0,1) );
  BOOST_CHECK_EQUAL(  14, sm2a(1,0) );
  BOOST_CHECK_EQUAL(  20, sm2a(1,1) );

  MatrixSymS2 sm2b = m1*Transpose(m1);

  BOOST_CHECK_EQUAL(   5, sm2b(0,0) );
  BOOST_CHECK_EQUAL(  11, sm2b(0,1) );
  BOOST_CHECK_EQUAL(  11, sm2b(1,0) );
  BOOST_CHECK_EQUAL(  25, sm2b(1,1) );


  VectorS2 v1 = {2,3};

  MatrixSymS2 sm3a = Transpose(m1)*diag(v1)*m1;

  BOOST_CHECK_EQUAL(  29, sm3a(0,0) );
  BOOST_CHECK_EQUAL(  40, sm3a(0,1) );
  BOOST_CHECK_EQUAL(  40, sm3a(1,0) );
  BOOST_CHECK_EQUAL(  56, sm3a(1,1) );

  sm3a += Transpose(m1)*diag(v1)*m1;

  BOOST_CHECK_EQUAL( 2*29, sm3a(0,0) );
  BOOST_CHECK_EQUAL( 2*40, sm3a(0,1) );
  BOOST_CHECK_EQUAL( 2*40, sm3a(1,0) );
  BOOST_CHECK_EQUAL( 2*56, sm3a(1,1) );

  sm3a -= Transpose(m1)*diag(v1)*m1;

  BOOST_CHECK_EQUAL(  29, sm3a(0,0) );
  BOOST_CHECK_EQUAL(  40, sm3a(0,1) );
  BOOST_CHECK_EQUAL(  40, sm3a(1,0) );
  BOOST_CHECK_EQUAL(  56, sm3a(1,1) );


  MatrixSymS2 sm3b = m1*diag(v1)*Transpose(m1);

  BOOST_CHECK_EQUAL(  14, sm3b(0,0) );
  BOOST_CHECK_EQUAL(  30, sm3b(0,1) );
  BOOST_CHECK_EQUAL(  30, sm3b(1,0) );
  BOOST_CHECK_EQUAL(  66, sm3b(1,1) );

  sm3b += m1*diag(v1)*Transpose(m1);

  BOOST_CHECK_EQUAL( 2*14, sm3b(0,0) );
  BOOST_CHECK_EQUAL( 2*30, sm3b(0,1) );
  BOOST_CHECK_EQUAL( 2*30, sm3b(1,0) );
  BOOST_CHECK_EQUAL( 2*66, sm3b(1,1) );

  sm3b -= m1*diag(v1)*Transpose(m1);

  BOOST_CHECK_EQUAL(  14, sm3b(0,0) );
  BOOST_CHECK_EQUAL(  30, sm3b(0,1) );
  BOOST_CHECK_EQUAL(  30, sm3b(1,0) );
  BOOST_CHECK_EQUAL(  66, sm3b(1,1) );


  MatrixSymS2 sm4 = {{2},{4,3}};

  MatrixSymS2 sm5a = Transpose(m1)*sm4*m1;
  BOOST_CHECK_EQUAL(  53, sm5a(0,0) );
  BOOST_CHECK_EQUAL(  80, sm5a(0,1) );
  BOOST_CHECK_EQUAL(  80, sm5a(1,0) );
  BOOST_CHECK_EQUAL( 120, sm5a(1,1) );

  sm5a += Transpose(m1)*sm4*m1;
  BOOST_CHECK_EQUAL(  2*53, sm5a(0,0) );
  BOOST_CHECK_EQUAL(  2*80, sm5a(0,1) );
  BOOST_CHECK_EQUAL(  2*80, sm5a(1,0) );
  BOOST_CHECK_EQUAL( 2*120, sm5a(1,1) );

  sm5a -= Transpose(m1)*sm4*m1;
  BOOST_CHECK_EQUAL(  53, sm5a(0,0) );
  BOOST_CHECK_EQUAL(  80, sm5a(0,1) );
  BOOST_CHECK_EQUAL(  80, sm5a(1,0) );
  BOOST_CHECK_EQUAL( 120, sm5a(1,1) );


  MatrixSymS2 sm5b = m1*sm4*Transpose(m1);
  BOOST_CHECK_EQUAL(  30, sm5b(0,0) );
  BOOST_CHECK_EQUAL(  70, sm5b(0,1) );
  BOOST_CHECK_EQUAL(  70, sm5b(1,0) );
  BOOST_CHECK_EQUAL( 162, sm5b(1,1) );

  sm5b += m1*sm4*Transpose(m1);
  BOOST_CHECK_EQUAL(  2*30, sm5b(0,0) );
  BOOST_CHECK_EQUAL(  2*70, sm5b(0,1) );
  BOOST_CHECK_EQUAL(  2*70, sm5b(1,0) );
  BOOST_CHECK_EQUAL( 2*162, sm5b(1,1) );

  sm5b -= m1*sm4*Transpose(m1);
  BOOST_CHECK_EQUAL(  30, sm5b(0,0) );
  BOOST_CHECK_EQUAL(  70, sm5b(0,1) );
  BOOST_CHECK_EQUAL(  70, sm5b(1,0) );
  BOOST_CHECK_EQUAL( 162, sm5b(1,1) );


  sm1a = {{2},
          {3, 4}};

  sm1b = {{5},
          {6, 7}};

  MatrixSymS2 sm6 = sm1a*sm1b*sm1a;

  BOOST_CHECK_EQUAL(  155, sm6(0,0) );
  BOOST_CHECK_EQUAL(  216, sm6(0,1) );
  BOOST_CHECK_EQUAL(  216, sm6(1,0) );
  BOOST_CHECK_EQUAL(  301, sm6(1,1) );

  sm6 += sm1a*sm1b*sm1a;

  BOOST_CHECK_EQUAL( 2*155, sm6(0,0) );
  BOOST_CHECK_EQUAL( 2*216, sm6(0,1) );
  BOOST_CHECK_EQUAL( 2*216, sm6(1,0) );
  BOOST_CHECK_EQUAL( 2*301, sm6(1,1) );

  sm6 -= sm1a*sm1b*sm1a;

  BOOST_CHECK_EQUAL(  155, sm6(0,0) );
  BOOST_CHECK_EQUAL(  216, sm6(0,1) );
  BOOST_CHECK_EQUAL(  216, sm6(1,0) );
  BOOST_CHECK_EQUAL(  301, sm6(1,1) );

  MatrixSymS2 sm7 = sm1a*sm1a;

  BOOST_CHECK_EQUAL( 13, sm7(0,0) );
  BOOST_CHECK_EQUAL( 18, sm7(0,1) );
  BOOST_CHECK_EQUAL( 18, sm7(1,0) );
  BOOST_CHECK_EQUAL( 25, sm7(1,1) );

  sm7 += sm1a*sm1a;

  BOOST_CHECK_EQUAL( 2*13, sm7(0,0) );
  BOOST_CHECK_EQUAL( 2*18, sm7(0,1) );
  BOOST_CHECK_EQUAL( 2*18, sm7(1,0) );
  BOOST_CHECK_EQUAL( 2*25, sm7(1,1) );

  sm7 -= sm1a*sm1a;

  BOOST_CHECK_EQUAL( 13, sm7(0,0) );
  BOOST_CHECK_EQUAL( 18, sm7(0,1) );
  BOOST_CHECK_EQUAL( 18, sm7(1,0) );
  BOOST_CHECK_EQUAL( 25, sm7(1,1) );


  //Make sure you can't construct a potentially non-symmetric matrix
  BOOST_CHECK_THROW( sm1a  = Transpose(m1)*m2, AssertionException );
  BOOST_CHECK_THROW( sm1a += Transpose(m1)*m2, AssertionException );
  BOOST_CHECK_THROW( sm1a -= Transpose(m1)*m2, AssertionException );

  BOOST_CHECK_THROW( sm1b  = Transpose(m1)*diag(v1)*m2, AssertionException );
  BOOST_CHECK_THROW( sm1b += Transpose(m1)*diag(v1)*m2, AssertionException );
  BOOST_CHECK_THROW( sm1b -= Transpose(m1)*diag(v1)*m2, AssertionException );

  BOOST_CHECK_THROW( sm1b  = Transpose(m1)*sm4*m2, AssertionException );
  BOOST_CHECK_THROW( sm1b += Transpose(m1)*sm4*m2, AssertionException );
  BOOST_CHECK_THROW( sm1b -= Transpose(m1)*sm4*m2, AssertionException );

  BOOST_CHECK_THROW( sm1b  = sm1a*sm1b*sm2a, AssertionException );
  BOOST_CHECK_THROW( sm1b += sm1a*sm1b*sm2a, AssertionException );
  BOOST_CHECK_THROW( sm1b -= sm1a*sm1b*sm2a, AssertionException );

  BOOST_CHECK_THROW( sm1b  = sm1a*sm1b, AssertionException );
  BOOST_CHECK_THROW( sm1b += sm1a*sm1b, AssertionException );
  BOOST_CHECK_THROW( sm1b -= sm1a*sm1b, AssertionException );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( ATransposeA_rectangle_ctor )
{

  const Int data[6] = {2,3,4,
                       5,6,7};
  MatrixS<2,3,Int> ma(data,6);

  MatrixSymS<2,Int> sm1a;
  sm1a = ma*Transpose(ma);
  BOOST_CHECK_EQUAL(  29, sm1a(0,0) );
  BOOST_CHECK_EQUAL(  56, sm1a(0,1) );
  BOOST_CHECK_EQUAL(  56, sm1a(1,0) );
  BOOST_CHECK_EQUAL( 110, sm1a(1,1) );

  sm1a += ma*Transpose(ma);

  BOOST_CHECK_EQUAL(  2*29, sm1a(0,0) );
  BOOST_CHECK_EQUAL(  2*56, sm1a(0,1) );
  BOOST_CHECK_EQUAL(  2*56, sm1a(1,0) );
  BOOST_CHECK_EQUAL( 2*110, sm1a(1,1) );

  sm1a -= ma*Transpose(ma);

  BOOST_CHECK_EQUAL(  29, sm1a(0,0) );
  BOOST_CHECK_EQUAL(  56, sm1a(0,1) );
  BOOST_CHECK_EQUAL(  56, sm1a(1,0) );
  BOOST_CHECK_EQUAL( 110, sm1a(1,1) );

  MatrixSymS<3,Int> sm1b;
  sm1b = Transpose(ma)*ma;

  BOOST_CHECK_EQUAL( 29, sm1b(0,0) );
  BOOST_CHECK_EQUAL( 36, sm1b(0,1) );
  BOOST_CHECK_EQUAL( 43, sm1b(0,2) );
  BOOST_CHECK_EQUAL( 36, sm1b(1,0) );
  BOOST_CHECK_EQUAL( 45, sm1b(1,1) );
  BOOST_CHECK_EQUAL( 54, sm1b(1,2) );
  BOOST_CHECK_EQUAL( 43, sm1b(2,0) );
  BOOST_CHECK_EQUAL( 54, sm1b(2,1) );
  BOOST_CHECK_EQUAL( 65, sm1b(2,2) );

  sm1b += Transpose(ma)*ma;

  BOOST_CHECK_EQUAL( 2*29, sm1b(0,0) );
  BOOST_CHECK_EQUAL( 2*36, sm1b(0,1) );
  BOOST_CHECK_EQUAL( 2*43, sm1b(0,2) );
  BOOST_CHECK_EQUAL( 2*36, sm1b(1,0) );
  BOOST_CHECK_EQUAL( 2*45, sm1b(1,1) );
  BOOST_CHECK_EQUAL( 2*54, sm1b(1,2) );
  BOOST_CHECK_EQUAL( 2*43, sm1b(2,0) );
  BOOST_CHECK_EQUAL( 2*54, sm1b(2,1) );
  BOOST_CHECK_EQUAL( 2*65, sm1b(2,2) );

  sm1b -= Transpose(ma)*ma;

  BOOST_CHECK_EQUAL( 29, sm1b(0,0) );
  BOOST_CHECK_EQUAL( 36, sm1b(0,1) );
  BOOST_CHECK_EQUAL( 43, sm1b(0,2) );
  BOOST_CHECK_EQUAL( 36, sm1b(1,0) );
  BOOST_CHECK_EQUAL( 45, sm1b(1,1) );
  BOOST_CHECK_EQUAL( 54, sm1b(1,2) );
  BOOST_CHECK_EQUAL( 43, sm1b(2,0) );
  BOOST_CHECK_EQUAL( 54, sm1b(2,1) );
  BOOST_CHECK_EQUAL( 65, sm1b(2,2) );

  const Int data2[6] = {7,6,
                        5,4,
                        3,2};
  MatrixS<3,2,Int> mb(data2,6);
  MatrixSymS<2,Int> sm2a;

  sm2a = Transpose(mb)*mb;

  BOOST_CHECK_EQUAL( 83, sm2a(0,0) );
  BOOST_CHECK_EQUAL( 68, sm2a(0,1) );
  BOOST_CHECK_EQUAL( 68, sm2a(1,0) );
  BOOST_CHECK_EQUAL( 56, sm2a(1,1) );

  sm2a += Transpose(mb)*mb;

  BOOST_CHECK_EQUAL( 2*83, sm2a(0,0) );
  BOOST_CHECK_EQUAL( 2*68, sm2a(0,1) );
  BOOST_CHECK_EQUAL( 2*68, sm2a(1,0) );
  BOOST_CHECK_EQUAL( 2*56, sm2a(1,1) );

  sm2a -= Transpose(mb)*mb;

  BOOST_CHECK_EQUAL( 83, sm2a(0,0) );
  BOOST_CHECK_EQUAL( 68, sm2a(0,1) );
  BOOST_CHECK_EQUAL( 68, sm2a(1,0) );
  BOOST_CHECK_EQUAL( 56, sm2a(1,1) );

  MatrixSymS<3,Int> sm2b;

  sm2b = mb*Transpose(mb);

  BOOST_CHECK_EQUAL( 85, sm2b(0,0) );
  BOOST_CHECK_EQUAL( 59, sm2b(0,1) );
  BOOST_CHECK_EQUAL( 33, sm2b(0,2) );
  BOOST_CHECK_EQUAL( 59, sm2b(1,0) );
  BOOST_CHECK_EQUAL( 41, sm2b(1,1) );
  BOOST_CHECK_EQUAL( 23, sm2b(1,2) );
  BOOST_CHECK_EQUAL( 33, sm2b(2,0) );
  BOOST_CHECK_EQUAL( 23, sm2b(2,1) );
  BOOST_CHECK_EQUAL( 13, sm2b(2,2) );

  sm2b += mb*Transpose(mb);

  BOOST_CHECK_EQUAL( 2*85, sm2b(0,0) );
  BOOST_CHECK_EQUAL( 2*59, sm2b(0,1) );
  BOOST_CHECK_EQUAL( 2*33, sm2b(0,2) );
  BOOST_CHECK_EQUAL( 2*59, sm2b(1,0) );
  BOOST_CHECK_EQUAL( 2*41, sm2b(1,1) );
  BOOST_CHECK_EQUAL( 2*23, sm2b(1,2) );
  BOOST_CHECK_EQUAL( 2*33, sm2b(2,0) );
  BOOST_CHECK_EQUAL( 2*23, sm2b(2,1) );
  BOOST_CHECK_EQUAL( 2*13, sm2b(2,2) );

  sm2b -= mb*Transpose(mb);

  BOOST_CHECK_EQUAL( 85, sm2b(0,0) );
  BOOST_CHECK_EQUAL( 59, sm2b(0,1) );
  BOOST_CHECK_EQUAL( 33, sm2b(0,2) );
  BOOST_CHECK_EQUAL( 59, sm2b(1,0) );
  BOOST_CHECK_EQUAL( 41, sm2b(1,1) );
  BOOST_CHECK_EQUAL( 23, sm2b(1,2) );
  BOOST_CHECK_EQUAL( 33, sm2b(2,0) );
  BOOST_CHECK_EQUAL( 23, sm2b(2,1) );
  BOOST_CHECK_EQUAL( 13, sm2b(2,2) );

  VectorS2 v1 = {1,1};
  VectorS<3,Int> v2 = {1,1,1};
  MatrixSymS2 sm3a;

  sm3a = ma*diag(v2)*Transpose(ma);

  BOOST_CHECK_EQUAL(  29, sm3a(0,0) );
  BOOST_CHECK_EQUAL(  56, sm3a(0,1) );
  BOOST_CHECK_EQUAL(  56, sm3a(1,0) );
  BOOST_CHECK_EQUAL( 110, sm3a(1,1) );

  sm3a += ma*diag(v2)*Transpose(ma);

  BOOST_CHECK_EQUAL(  2*29, sm3a(0,0) );
  BOOST_CHECK_EQUAL(  2*56, sm3a(0,1) );
  BOOST_CHECK_EQUAL(  2*56, sm3a(1,0) );
  BOOST_CHECK_EQUAL( 2*110, sm3a(1,1) );

  sm3a -= ma*diag(v2)*Transpose(ma);

  BOOST_CHECK_EQUAL(  29, sm3a(0,0) );
  BOOST_CHECK_EQUAL(  56, sm3a(0,1) );
  BOOST_CHECK_EQUAL(  56, sm3a(1,0) );
  BOOST_CHECK_EQUAL( 110, sm3a(1,1) );

  MatrixSymS<3,Int> sm3b;
  sm3b = Transpose(ma)*diag(v1)*ma;

  BOOST_CHECK_EQUAL( 29, sm3b(0,0) );
  BOOST_CHECK_EQUAL( 36, sm3b(0,1) );
  BOOST_CHECK_EQUAL( 43, sm3b(0,2) );
  BOOST_CHECK_EQUAL( 36, sm3b(1,0) );
  BOOST_CHECK_EQUAL( 45, sm3b(1,1) );
  BOOST_CHECK_EQUAL( 54, sm3b(1,2) );
  BOOST_CHECK_EQUAL( 43, sm3b(2,0) );
  BOOST_CHECK_EQUAL( 54, sm3b(2,1) );
  BOOST_CHECK_EQUAL( 65, sm3b(2,2) );

  sm3b += Transpose(ma)*diag(v1)*ma;

  BOOST_CHECK_EQUAL( 2*29, sm3b(0,0) );
  BOOST_CHECK_EQUAL( 2*36, sm3b(0,1) );
  BOOST_CHECK_EQUAL( 2*43, sm3b(0,2) );
  BOOST_CHECK_EQUAL( 2*36, sm3b(1,0) );
  BOOST_CHECK_EQUAL( 2*45, sm3b(1,1) );
  BOOST_CHECK_EQUAL( 2*54, sm3b(1,2) );
  BOOST_CHECK_EQUAL( 2*43, sm3b(2,0) );
  BOOST_CHECK_EQUAL( 2*54, sm3b(2,1) );
  BOOST_CHECK_EQUAL( 2*65, sm3b(2,2) );

  sm3b -= Transpose(ma)*diag(v1)*ma;

  BOOST_CHECK_EQUAL( 29, sm3b(0,0) );
  BOOST_CHECK_EQUAL( 36, sm3b(0,1) );
  BOOST_CHECK_EQUAL( 43, sm3b(0,2) );
  BOOST_CHECK_EQUAL( 36, sm3b(1,0) );
  BOOST_CHECK_EQUAL( 45, sm3b(1,1) );
  BOOST_CHECK_EQUAL( 54, sm3b(1,2) );
  BOOST_CHECK_EQUAL( 43, sm3b(2,0) );
  BOOST_CHECK_EQUAL( 54, sm3b(2,1) );
  BOOST_CHECK_EQUAL( 65, sm3b(2,2) );

  MatrixSymS<2,Int> sm4a;

  sm4a = Transpose(mb)*diag(v2)*mb;

  BOOST_CHECK_EQUAL( 83, sm2a(0,0) );
  BOOST_CHECK_EQUAL( 68, sm2a(0,1) );
  BOOST_CHECK_EQUAL( 68, sm2a(1,0) );
  BOOST_CHECK_EQUAL( 56, sm2a(1,1) );

  sm4a += Transpose(mb)*diag(v2)*mb;

  BOOST_CHECK_EQUAL( 2*83, sm4a(0,0) );
  BOOST_CHECK_EQUAL( 2*68, sm4a(0,1) );
  BOOST_CHECK_EQUAL( 2*68, sm4a(1,0) );
  BOOST_CHECK_EQUAL( 2*56, sm4a(1,1) );

  sm4a -= Transpose(mb)*diag(v2)*mb;

  BOOST_CHECK_EQUAL( 83, sm2a(0,0) );
  BOOST_CHECK_EQUAL( 68, sm2a(0,1) );
  BOOST_CHECK_EQUAL( 68, sm2a(1,0) );
  BOOST_CHECK_EQUAL( 56, sm2a(1,1) );

  MatrixSymS<3,Int> sm4b;

  sm4b = mb*diag(v1)*Transpose(mb);

  BOOST_CHECK_EQUAL( 85, sm4b(0,0) );
  BOOST_CHECK_EQUAL( 59, sm4b(0,1) );
  BOOST_CHECK_EQUAL( 33, sm4b(0,2) );
  BOOST_CHECK_EQUAL( 59, sm4b(1,0) );
  BOOST_CHECK_EQUAL( 41, sm4b(1,1) );
  BOOST_CHECK_EQUAL( 23, sm4b(1,2) );
  BOOST_CHECK_EQUAL( 33, sm4b(2,0) );
  BOOST_CHECK_EQUAL( 23, sm4b(2,1) );
  BOOST_CHECK_EQUAL( 13, sm4b(2,2) );

  sm4b += mb*diag(v1)*Transpose(mb);

  BOOST_CHECK_EQUAL( 2*85, sm4b(0,0) );
  BOOST_CHECK_EQUAL( 2*59, sm4b(0,1) );
  BOOST_CHECK_EQUAL( 2*33, sm4b(0,2) );
  BOOST_CHECK_EQUAL( 2*59, sm4b(1,0) );
  BOOST_CHECK_EQUAL( 2*41, sm4b(1,1) );
  BOOST_CHECK_EQUAL( 2*23, sm4b(1,2) );
  BOOST_CHECK_EQUAL( 2*33, sm4b(2,0) );
  BOOST_CHECK_EQUAL( 2*23, sm4b(2,1) );
  BOOST_CHECK_EQUAL( 2*13, sm4b(2,2) );

  sm4b -= mb*diag(v1)*Transpose(mb);

  BOOST_CHECK_EQUAL( 85, sm4b(0,0) );
  BOOST_CHECK_EQUAL( 59, sm4b(0,1) );
  BOOST_CHECK_EQUAL( 33, sm4b(0,2) );
  BOOST_CHECK_EQUAL( 59, sm4b(1,0) );
  BOOST_CHECK_EQUAL( 41, sm4b(1,1) );
  BOOST_CHECK_EQUAL( 23, sm4b(1,2) );
  BOOST_CHECK_EQUAL( 33, sm4b(2,0) );
  BOOST_CHECK_EQUAL( 23, sm4b(2,1) );
  BOOST_CHECK_EQUAL( 13, sm4b(2,2) );

  MatrixSymS2 sma = {{2},{4,3}};

  MatrixSymS<3,Int> sm5a;

  sm5a = Transpose(ma)*sma*ma;

  BOOST_CHECK_EQUAL( 163, sm5a(0,0) );
  BOOST_CHECK_EQUAL( 210, sm5a(0,1) );
  BOOST_CHECK_EQUAL( 257, sm5a(0,2) );
  BOOST_CHECK_EQUAL( 210, sm5a(1,0) );
  BOOST_CHECK_EQUAL( 270, sm5a(1,1) );
  BOOST_CHECK_EQUAL( 330, sm5a(1,2) );
  BOOST_CHECK_EQUAL( 257, sm5a(2,0) );
  BOOST_CHECK_EQUAL( 330, sm5a(2,1) );
  BOOST_CHECK_EQUAL( 403, sm5a(2,2) );

  sm5a += Transpose(ma)*sma*ma;

  BOOST_CHECK_EQUAL( 2*163, sm5a(0,0) );
  BOOST_CHECK_EQUAL( 2*210, sm5a(0,1) );
  BOOST_CHECK_EQUAL( 2*257, sm5a(0,2) );
  BOOST_CHECK_EQUAL( 2*210, sm5a(1,0) );
  BOOST_CHECK_EQUAL( 2*270, sm5a(1,1) );
  BOOST_CHECK_EQUAL( 2*330, sm5a(1,2) );
  BOOST_CHECK_EQUAL( 2*257, sm5a(2,0) );
  BOOST_CHECK_EQUAL( 2*330, sm5a(2,1) );
  BOOST_CHECK_EQUAL( 2*403, sm5a(2,2) );

  sm5a -= Transpose(ma)*sma*ma;

  BOOST_CHECK_EQUAL( 163, sm5a(0,0) );
  BOOST_CHECK_EQUAL( 210, sm5a(0,1) );
  BOOST_CHECK_EQUAL( 257, sm5a(0,2) );
  BOOST_CHECK_EQUAL( 210, sm5a(1,0) );
  BOOST_CHECK_EQUAL( 270, sm5a(1,1) );
  BOOST_CHECK_EQUAL( 330, sm5a(1,2) );
  BOOST_CHECK_EQUAL( 257, sm5a(2,0) );
  BOOST_CHECK_EQUAL( 330, sm5a(2,1) );
  BOOST_CHECK_EQUAL( 403, sm5a(2,2) );

  MatrixSymS<3,Int> smb = {{2},
                           {4,3},
                           {4,5,8}};

  MatrixSymS<2,Int> sm5b;

  sm5b = ma*smb*Transpose(ma);

  BOOST_CHECK_EQUAL(  395, sm5b(0,0) );
  BOOST_CHECK_EQUAL(  767, sm5b(0,1) );
  BOOST_CHECK_EQUAL(  767, sm5b(1,0) );
  BOOST_CHECK_EQUAL( 1490, sm5b(1,1) );

  sm5b += ma*smb*Transpose(ma);

  BOOST_CHECK_EQUAL(  2*395, sm5b(0,0) );
  BOOST_CHECK_EQUAL(  2*767, sm5b(0,1) );
  BOOST_CHECK_EQUAL(  2*767, sm5b(1,0) );
  BOOST_CHECK_EQUAL( 2*1490, sm5b(1,1) );

  sm5b -= ma*smb*Transpose(ma);

  BOOST_CHECK_EQUAL(  395, sm5b(0,0) );
  BOOST_CHECK_EQUAL(  767, sm5b(0,1) );
  BOOST_CHECK_EQUAL(  767, sm5b(1,0) );
  BOOST_CHECK_EQUAL( 1490, sm5b(1,1) );

  MatrixSymS<3,Int> sm6a;

  sm6a = mb*sma*Transpose(mb);

  BOOST_CHECK_EQUAL( 542, sm6a(0,0) );
  BOOST_CHECK_EQUAL( 374, sm6a(0,1) );
  BOOST_CHECK_EQUAL( 206, sm6a(0,2) );
  BOOST_CHECK_EQUAL( 374, sm6a(1,0) );
  BOOST_CHECK_EQUAL( 258, sm6a(1,1) );
  BOOST_CHECK_EQUAL( 142, sm6a(1,2) );
  BOOST_CHECK_EQUAL( 206, sm6a(2,0) );
  BOOST_CHECK_EQUAL( 142, sm6a(2,1) );
  BOOST_CHECK_EQUAL(  78, sm6a(2,2) );

  sm6a += mb*sma*Transpose(mb);

  BOOST_CHECK_EQUAL( 2*542, sm6a(0,0) );
  BOOST_CHECK_EQUAL( 2*374, sm6a(0,1) );
  BOOST_CHECK_EQUAL( 2*206, sm6a(0,2) );
  BOOST_CHECK_EQUAL( 2*374, sm6a(1,0) );
  BOOST_CHECK_EQUAL( 2*258, sm6a(1,1) );
  BOOST_CHECK_EQUAL( 2*142, sm6a(1,2) );
  BOOST_CHECK_EQUAL( 2*206, sm6a(2,0) );
  BOOST_CHECK_EQUAL( 2*142, sm6a(2,1) );
  BOOST_CHECK_EQUAL(  2*78, sm6a(2,2) );

  sm6a -= mb*sma*Transpose(mb);

  BOOST_CHECK_EQUAL( 542, sm6a(0,0) );
  BOOST_CHECK_EQUAL( 374, sm6a(0,1) );
  BOOST_CHECK_EQUAL( 206, sm6a(0,2) );
  BOOST_CHECK_EQUAL( 374, sm6a(1,0) );
  BOOST_CHECK_EQUAL( 258, sm6a(1,1) );
  BOOST_CHECK_EQUAL( 142, sm6a(1,2) );
  BOOST_CHECK_EQUAL( 206, sm6a(2,0) );
  BOOST_CHECK_EQUAL( 142, sm6a(2,1) );
  BOOST_CHECK_EQUAL(  78, sm6a(2,2) );

  MatrixSymS<2,Int> sm6b;

  sm6b = Transpose(mb)*smb*mb;

  BOOST_CHECK_EQUAL( 843, sm6b(0,0) );
  BOOST_CHECK_EQUAL( 662, sm6b(0,1) );
  BOOST_CHECK_EQUAL( 662, sm6b(1,0) );
  BOOST_CHECK_EQUAL( 520, sm6b(1,1) );

  sm6b += Transpose(mb)*smb*mb;

  BOOST_CHECK_EQUAL( 2*843, sm6b(0,0) );
  BOOST_CHECK_EQUAL( 2*662, sm6b(0,1) );
  BOOST_CHECK_EQUAL( 2*662, sm6b(1,0) );
  BOOST_CHECK_EQUAL( 2*520, sm6b(1,1) );

  sm6b -= Transpose(mb)*smb*mb;

  BOOST_CHECK_EQUAL( 843, sm6b(0,0) );
  BOOST_CHECK_EQUAL( 662, sm6b(0,1) );
  BOOST_CHECK_EQUAL( 662, sm6b(1,0) );
  BOOST_CHECK_EQUAL( 520, sm6b(1,1) );

  MatrixS<3,2,Int> mc(data,6);
  MatrixS<2,3,Int> md(data,6);

  //Make sure you can't construct a potentially non-symmetric matrix
  BOOST_CHECK_THROW( sm1a  = Transpose(mb)*mc;, AssertionException );
  BOOST_CHECK_THROW( sm1a += Transpose(mb)*mc;, AssertionException );
  BOOST_CHECK_THROW( sm1a -= Transpose(mb)*mc;, AssertionException );

  BOOST_CHECK_THROW( sm1b  = Transpose(md)*diag(v1)*ma;, AssertionException );
  BOOST_CHECK_THROW( sm1b += Transpose(md)*diag(v1)*ma;, AssertionException );
  BOOST_CHECK_THROW( sm1b -= Transpose(md)*diag(v1)*ma;, AssertionException );

  BOOST_CHECK_THROW( sm1b  = Transpose(md)*sm6b*ma;, AssertionException );
  BOOST_CHECK_THROW( sm1b += Transpose(md)*sm6b*ma;, AssertionException );
  BOOST_CHECK_THROW( sm1b -= Transpose(md)*sm6b*ma;, AssertionException );
}


//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( MatrixSymS_Transpose_xTMx_test )
{
  Int Mdata[] = {1,
                 3, 4};

  Int xdata[] = {3, 4};

  MatrixSymS2 M(Mdata, 3);
  VectorS2 x(xdata, 2);
  Int xTMx;

  xTMx = Transpose(x)*M*x;

  BOOST_CHECK_EQUAL( xTMx , 145 );
}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( matrix_ops2 )
{
  MatrixSymS2 m1 = {{1},
                    {2,4}};
  MatrixSymS2 m2(m1), m3, m4, m5;

  // size
  BOOST_CHECK( m1.M == 2 );
  BOOST_CHECK( m1.N == 2 );
  BOOST_CHECK( m2.M == 2 );
  BOOST_CHECK( m2.N == 2 );

  // ctors
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );

  // assignment
  m3 = m1;
  m4 = 5;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 5,5,5,5 ) );

  m2 = m3 = 3;
  BOOST_CHECK( chkMatrixSymS22( m2, 3,3,3,3 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 3,3,3,3 ) );

  // unary
  m2 = +m1;
  m3 = -m1;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );

  // binary accumulation
  m4 = m1;
  m4 *= 5;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 5,10,10,20 ) );

  m2 = 5;
  m3 = 5;
  m2 += m1;
  m3 -= m1;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 6,7,7,9 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 4,3,3,1 ) );

  // binary operators
//  m2 = m1 + 3;
//  m3 = m1 - 3;
  m4 = m1 * 3;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
//  BOOST_CHECK( chkMatrixSymS22( m2, 4,5,6,7 ) );
//  BOOST_CHECK( chkMatrixSymS22( m3, -2,-1,0,1 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 3,6,6,12 ) );

//  m2 = 3 + m1;
//  m3 = 3 - m1;
  m4 = 3 * m1;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
//  BOOST_CHECK( chkMatrixSymS22( m2, 4,5,6,7 ) );
//  BOOST_CHECK( chkMatrixSymS22( m3, 2,1,0,-1 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 3,6,6,12 ) );

  m2 = 3;
  m3 = m1 + m2;
  m4 = m1 - m2;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 3,3,3,3 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 4,5,5,7 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, -2,-1,-1,1 ) );

  // arithmetic combinations

  m2 = m1;
  m3 = m1 + m2;
  m4 = m1 + m2 + m3;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 2,4,4,8 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 4,8,8,16 ) );

  m2 += m1;
  m3 += m1 + m2;
  m4 += m1 + m2 + m3;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 2,4,4,8 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 5,10,10,20 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 12,24,24,48 ) );

  m3 = m1 - m2;
  m4 = m1 - m2 - m3;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 2,4,4,8 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 0,0,0,0 ) );

  m2 -= m1;
  m3 -= m1 - m2;
  m4 -= m1 - m2 - m3;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );

  m3 = m1 - m2;
  m4 = m1 + m2 - m3;
  m5 = m1 - m2 + m3;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 0,0,0,0 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );
  BOOST_CHECK( chkMatrixSymS22( m5, 0,0,0,0 ) );

  m5 = (m1 + m2) + (m3 + m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = (m1 + m2) + (m3 - m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 0,0,0,0 ) );
  m5 = (m1 + m2) - (m3 + m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 0,0,0,0 ) );
  m5 = (m1 + m2) - (m3 - m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = (m1 - m2) + (m3 + m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 2,4,4,8 ) );
  m5 = (m1 - m2) + (m3 - m4);
  BOOST_CHECK( chkMatrixSymS22( m5, -2,-4,-4,-8 ) );
  m5 = (m1 - m2) - (m3 + m4);
  BOOST_CHECK( chkMatrixSymS22( m5, -2,-4,-4,-8 ) );
  m5 = (m1 - m2) - (m3 - m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 2,4,4,8 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 0,0,0,0 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 += (m1 + m2) + (m3 + m4);
  m5 += (m1 + m2) + (m3 - m4);
  m5 += (m1 + m2) - (m3 + m4);
  m5 += (m1 + m2) - (m3 - m4);
  m5 += (m1 - m2) + (m3 + m4);
  m5 += (m1 - m2) + (m3 - m4);
  m5 += (m1 - m2) - (m3 + m4);
  m5 += (m1 - m2) - (m3 - m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 10,20,20,40 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 0,0,0,0 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 -= (m1 + m2) + (m3 + m4);
  m5 -= (m1 + m2) + (m3 - m4);
  m5 -= (m1 + m2) - (m3 + m4);
  m5 -= (m1 + m2) - (m3 - m4);
  m5 -= (m1 - m2) + (m3 + m4);
  m5 -= (m1 - m2) + (m3 - m4);
  m5 -= (m1 - m2) - (m3 + m4);
  m5 -= (m1 - m2) - (m3 - m4);
  BOOST_CHECK( chkMatrixSymS22( m5, 2,4,4,8 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 0,0,0,0 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m2 = 1*m1;
  m3 = m2*2;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 2,4,4,8 ) );

  m2 += 1*m1;
  m3 += m2*2;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 2,4,4,8 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 6,12,12,24 ) );

  m2 -= 1*m1;
  m3 -= m2*2;
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 4,8,8,16 ) );



  m3 = (2*m2)*2;
  BOOST_CHECK( chkMatrixSymS22( m3, 4,8,8,16 ) );

  m3 += (2*m2)*2;
  BOOST_CHECK( chkMatrixSymS22( m3, 2*4,2*8,2*8,2*16 ) );

  m3 = 2*(m2*2);
  BOOST_CHECK( chkMatrixSymS22( m3, 4,8,8,16 ) );

  m3 += 2*(m2*2);
  BOOST_CHECK( chkMatrixSymS22( m3, 2*4,2*8,2*8,2*16 ) );

  m3 = (4*m2)/2;
  BOOST_CHECK( chkMatrixSymS22( m3, 2,4,4,8 ) );

  m3 += (4*m2)/2;
  BOOST_CHECK( chkMatrixSymS22( m3, 2*2,2*4,2*4,2*8 ) );



  m5 = 1*(m1 + m2) + (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 14,28,28,56 ) );
  m5 = 1*(m1 + m2) + (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 6,12,12,24 ) );
  m5 = 1*(m1 + m2) - (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -10,-20,-20,-40 ) );
  m5 = 1*(m1 + m2) - (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -2,-4,-4,-8 ) );
  m5 = 1*(m1 - m2) + (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 12,24,24,48 ) );
  m5 = 1*(m1 - m2) + (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = 1*(m1 - m2) - (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -12,-24,-24,-48 ) );
  m5 = 1*(m1 - m2) - (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -4,-8,-8,-16 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 4,8,8,16 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 += 1*(m1 + m2) + (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 10,20,20,40 ) );
  m5 += 1*(m1 + m2) + (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 16,32,32,64 ) );
  m5 += 1*(m1 + m2) - (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 6,12,12,24 ) );
  m5 += 1*(m1 + m2) - (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 += 1*(m1 - m2) + (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 16,32,32,64 ) );
  m5 += 1*(m1 - m2) + (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 20,40,40,80 ) );
  m5 += 1*(m1 - m2) - (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 8,16,16,32 ) );
  m5 += 1*(m1 - m2) - (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 4,8,8,16 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 -= 1*(m1 + m2) + (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -10,-20,-20,-40 ) );
  m5 -= 1*(m1 + m2) + (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -16,-32,-32,-64 ) );
  m5 -= 1*(m1 + m2) - (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -6,-12,-12,-24 ) );
  m5 -= 1*(m1 + m2) - (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -4,-8,-8,-16 ) );
  m5 -= 1*(m1 - m2) + (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -16,-32,-32,-64 ) );
  m5 -= 1*(m1 - m2) + (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -20,-40,-40,-80 ) );
  m5 -= 1*(m1 - m2) - (m3 + m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -8,-16,-16,-32 ) );
  m5 -= 1*(m1 - m2) - (m3 - m4)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, -4,-8,-8,-16 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, 4,8,8,16 ) );
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 = 1*(m1 + m2)*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = 1*2*(m1 + m2);
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = (m1 + m2)*1*2;
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );

  //Need to divide by 1 because we are working with integers
  m5 = (m1 + m2)/1;
  BOOST_CHECK( chkMatrixSymS22( m5, 2,4,4,8 ) );

  m5 = +( 2*(m1 + m2) );
  BOOST_CHECK( chkMatrixSymS22( m5, 4,8,8,16 ) );

  m2 = +m1;
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  m3 = -m2;
  BOOST_CHECK( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );
  m4 = +(m1 + m2);
  BOOST_CHECK( chkMatrixSymS22( m4, 2,4,4,8 ) );
  m4 = +(m1 - m2);
  BOOST_CHECK( chkMatrixSymS22( m4, 0,0,0,0 ) );
  m4 = -(m1 + m2);
  BOOST_CHECK( chkMatrixSymS22( m4, -2,-4,-4,-8 ) );
  m4 = -(m1 - m2);
  BOOST_CHECK( chkMatrixSymS22( m4, 0,0,0,0 ) );
  m4 = +(m1 + m2) + m3;
  BOOST_CHECK( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -(m1 + m2) + m3;
  BOOST_CHECK( chkMatrixSymS22( m4, -3,-6,-6,-12 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 1,2,2,4 ) );
  BOOST_CHECK( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );

  m4 = +1*m1;
  BOOST_CHECK( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -1*m1;
  BOOST_CHECK( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  m4 = +m1*1;
  BOOST_CHECK( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -m1*1;
  BOOST_CHECK( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  m4 = +(1*m1);
  BOOST_CHECK( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -(1*m1);
  BOOST_CHECK( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  m4 = +(m1*1);
  BOOST_CHECK( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -(m1*1);
  BOOST_CHECK( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 1,2,2,4 ) );
}


//----------------------------------------------------------------------------//
// matrix-vector multiply
BOOST_AUTO_TEST_CASE( matrix_vector_multiply2 )
{
  const Int adata[2] = {1,2};
  const Int mdata[3] = {3,4,6};
  MatrixSymS2 m1(mdata, 3);
  MatrixSymS2 m2(m1), m3;
  VectorS2 a1(adata, 2);
  VectorS2 a2(a1), a3;

  BOOST_CHECK( chkVectorS2( a1, 1,2 ) );
  BOOST_CHECK( chkVectorS2( a2, 1,2 ) );

  BOOST_CHECK( chkMatrixSymS22( m1, 3,4,4,6 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 3,4,4,6 ) );

  a2 = m1*a1;
  BOOST_CHECK( chkVectorS2( a2, 11,16 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 3,4,4,6 ) );

  a2 = +(m1*a1);
  BOOST_CHECK( chkVectorS2( a2, 11,16 ) );

  a2 += m1*a1;
  BOOST_CHECK( chkVectorS2( a2, 2*11,2*16 ) );

  a2 -= m1*a1;
  BOOST_CHECK( chkVectorS2( a2, 11,16 ) );

  a2 += m1*a1 + m1*a1;
  BOOST_CHECK( chkVectorS2( a2, 3*11,3*16 ) );

  a2 += m1*a1 - m1*a1;
  BOOST_CHECK( chkVectorS2( a2, 3*11,3*16 ) );



  a2 = (2*m1)*a1;
  BOOST_CHECK( chkVectorS2( a2, 2*11,2*16 ) );

  a2 = +((2*m1)*a1);
  BOOST_CHECK( chkVectorS2( a2, 2*11,2*16 ) );

  a2 += (2*m1)*a1;
  BOOST_CHECK( chkVectorS2( a2, 4*11,4*16 ) );

  a2 = m1*(a1*2);
  BOOST_CHECK( chkVectorS2( a2, 2*11,2*16 ) );

  a2 = +(m1*(a1*2));
  BOOST_CHECK( chkVectorS2( a2, 2*11,2*16 ) );

  a2 += m1*(a1*2);
  BOOST_CHECK( chkVectorS2( a2, 4*11,4*16 ) );



  a2 = (m1 + m2)*a1;
  BOOST_CHECK( chkVectorS2( a2, 22,32 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 3,4,4,6 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 3,4,4,6 ) );

  a2 = +((m1 + m2)*a1);
  BOOST_CHECK( chkVectorS2( a2, 22,32 ) );

  a2 += (m1 + m2)*a1;
  BOOST_CHECK( chkVectorS2( a2, 2*22,2*32 ) );

  a2 -= (m1 + m2)*a1;
  BOOST_CHECK( chkVectorS2( a2, 22,32 ) );

//  a2 = a1*(m1 + m2);
//  BOOST_CHECK( chkVectorS2( a2, 26,32 ) );
//  BOOST_CHECK( chkMatrixSymS22( m1, 3,4,5,6 ) );
//  BOOST_CHECK( chkMatrixSymS22( m2, 3,4,5,6 ) );

  a2 = (m1 - m2)*a1;
  BOOST_CHECK( chkVectorS2( a2, 0,0 ) );
  BOOST_CHECK( chkMatrixSymS22( m1, 3,4,4,6 ) );
  BOOST_CHECK( chkMatrixSymS22( m2, 3,4,4,6 ) );

//  a2 = a1*(m1 - m2);
//  BOOST_CHECK( chkVectorS2( a2, 0,0 ) );
//  BOOST_CHECK( chkMatrixSymS22( m1, 3,4,5,6 ) );
//  BOOST_CHECK( chkMatrixSymS22( m2, 3,4,5,6 ) );

  a2 = a1;

  a3 = m1*(a1 + a2);
  BOOST_CHECK( chkVectorS2( a3, 2*11,2*16 ) );

  a3 = +(m1*(a1 + a2));
  BOOST_CHECK( chkVectorS2( a3, 2*11,2*16 ) );

  a3 += m1*(a1 + a2);
  BOOST_CHECK( chkVectorS2( a3, 4*11,4*16 ) );

  a3 -= m1*(a1 + a2);
  BOOST_CHECK( chkVectorS2( a3, 2*11,2*16 ) );


  a2 = a1;

  a3 = (m1 + m2)*(a1 + a2);
  BOOST_CHECK( chkVectorS2( a3, 44,64 ) );

  a3 += (m1 + m2)*(a1 + a2);
  BOOST_CHECK( chkVectorS2( a3, 2*44,2*64 ) );

  a3 -= (m1 + m2)*(a1 + a2);
  BOOST_CHECK( chkVectorS2( a3, 44,64 ) );

  a3 = +((m1 + m2)*(a1 + a2));
  BOOST_CHECK( chkVectorS2( a3, 44,64 ) );

//  a3 = (a1 + a2)*(m1 + m2);
//  BOOST_CHECK( chkVectorS2( a3, 52,64 ) );

  a3 = (m1 + m2)*(a1 - a2);
  BOOST_CHECK( chkVectorS2( a3, 0,0 ) );

//  a3 = (a1 - a2)*(m1 + m2);
//  BOOST_CHECK( chkVectorS2( a3, 0,0 ) );

  a3 = (m1 - m2)*(a1 + a2);
  BOOST_CHECK( chkVectorS2( a3, 0,0 ) );

//  a3 = (a1 + a2)*(m1 - m2);
//  BOOST_CHECK( chkVectorS2( a3, 0,0 ) );

  a3 = (m1 - m2)*(a1 - a2);
  BOOST_CHECK( chkVectorS2( a3, 0,0 ) );

//  a3 = (a1 - a2)*(m1 - m2);
//  BOOST_CHECK( chkVectorS2( a3, 0,0 ) );
//  BOOST_CHECK( chkVectorS2( a1, 1,2 ) );
//  BOOST_CHECK( chkVectorS2( a2, 1,2 ) );
//  BOOST_CHECK( chkMatrixSymS22( m1, 3,4,5,6 ) );
//  BOOST_CHECK( chkMatrixSymS22( m2, 3,4,5,6 ) );

}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( IO )
{
  //Set the 2nd argument to false to regenerate the pattern file
  output_test_stream output( "IO/DenseLinAlg/MatrixSymS_pattern.txt", true );

  MatrixSymS2 m = Identity();

  output << m << std::endl;
  BOOST_CHECK( output.match_pattern() );
  m.dump( 2, output );
  BOOST_CHECK( output.match_pattern() );
}

#endif

//############################################################################//
UT_TEST_SUITE_END(MatrixSymD_tests)
