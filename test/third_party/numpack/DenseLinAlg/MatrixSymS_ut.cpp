// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// MatrixS_btest
// testing of MatrixS<M,N,T> class

#include "unit_tester.hpp"

#include "numpack/dense/static/MatrixS.h"
//#include "numpack/dense/static/MatrixS_Sub.h"
#include "numpack/dense/static/MatrixSymS.h"
#include "numpack/dense/static/VectorS.h"

#include "chkMatrixD_btest.h"

#include <iostream>

using namespace numpack;
using namespace numpack::DLA;

UT_TEST_SUITE( MatrixSymD_tests )

typedef double Int;

/*typedef VectorS<2,Int> VectorS;
typedef MatrixS<2,2,Int> MatrixS;
typedef MatrixSymS<2,Int> MatrixSymS;*/

UT_TEST_CASE( test1 )
{
  const int n = 2;

  DLA::MatrixSymS<2,Real> m1;
  DLA::MatrixSymS<2,Real> m2;

  m1 = 1;
  m2 = 2;

  DLA::MatrixSymS<2,Real> m3;
  m3 =  m1+m2;
  m3.dump();

  DLA::MatrixS<2,2,Real> m4;

  m4 = m1*m2;
  m4.dump();
}
UT_TEST_CASE_END( test1 )

//----------------------------------------------------------------------------//
UT_TEST_CASE( matrix_ops1 )
{
  const Int data = 3.;
  MatrixSymS<1,double> m1(&data, 1);
  MatrixSymS<1,double> m2(m1);
  MatrixSymS<1,double> m3, m4, m5;

  m3 = 1;
  m4 = 2;
  m5 = 0;

  // size
  UT_ASSERT( m1.M == 1 );
  UT_ASSERT( m1.N == 1 );
  UT_ASSERT( m2.M == 1 );
  UT_ASSERT( m2.N == 1 );

  // ctors
  UT_ASSERT_EQUALS(  3, m1(0,0) );
  UT_ASSERT_EQUALS(  3, m2(0,0) );
  UT_ASSERT_EQUALS(  1, m3(0,0) );
  UT_ASSERT_EQUALS(  2, m4(0,0) );

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
  m4 = m1;
  m3 *= data;
  m4 /= data;
  UT_ASSERT_EQUALS(  9, m3(0,0) );
  UT_ASSERT_EQUALS(  1, m4(0,0) );

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
//  m2 = m1 + data;
//  m3 = m1 - data;
  m4 = m1 * data;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
//  UT_ASSERT_EQUALS(  6, m2(0,0) );
//  UT_ASSERT_EQUALS(  0, m3(0,0) );
  UT_ASSERT_EQUALS(  9, m4(0,0) );

//  m2 = data + m1;
//  m3 = data - m1;
  m4 = data * m1;
  UT_ASSERT_EQUALS(  3, m1(0,0) );
//  UT_ASSERT_EQUALS(  6, m2(0,0) );
//  UT_ASSERT_EQUALS(  0, m3(0,0) );
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

  m4 = {6};
  UT_ASSERT_EQUALS( 6, m4(0,0) );

  Int i0 = m1*m1;
  UT_ASSERT_EQUALS( 9, i0 );

  /*Int i1 = (m1+m1)*m1;
  UT_ASSERT_EQUALS( 18, i1 );

  Int i2 = m1*(m1+m1);
  UT_ASSERT_EQUALS( 18, i2 );

  Int i3 = (m1+m1)*(m1+m1);
  UT_ASSERT_EQUALS( 36, i3 );*/
}
UT_TEST_CASE_END( matrix_ops1 )

//----------------------------------------------------------------------------//
UT_TEST_CASE( identity2 )
{
  MatrixSymS<2,double> m1;
  m1 = Identity();

  UT_ASSERT_EQUALS(  1, m1(0,0) );
  UT_ASSERT_EQUALS(  0, m1(0,1) );
  UT_ASSERT_EQUALS(  0, m1(1,0) );
  UT_ASSERT_EQUALS(  1, m1(1,1) );

  MatrixSymS<2,double> m2;
  MatrixSymS<2,double> I;
  I = Identity();

  m2 = 2*I;

  UT_ASSERT_EQUALS(  2, m2(0,0) );
  UT_ASSERT_EQUALS(  0, m2(0,1) );
  UT_ASSERT_EQUALS(  0, m2(1,0) );
  UT_ASSERT_EQUALS(  2, m2(1,1) );

  MatrixSymS<2,double> m3;
  m3 = Identity();

  UT_ASSERT_EQUALS(  1, m3(0,0) );
  UT_ASSERT_EQUALS(  0, m3(0,1) );
  UT_ASSERT_EQUALS(  0, m3(1,0) );
  UT_ASSERT_EQUALS(  1, m3(1,1) );

}
UT_TEST_CASE_END(identity2)



//----------------------------------------------------------------------------//
UT_TEST_CASE( MatrixSymD_3x3_test )
{
  Int a = 1; Int b = 4; Int c = 7;
  Int d = 4; Int e = 5; Int f = 8;
  Int g = 7; Int h = 8; Int i = 9;

  Int Adata[] = {a,
                 d, e,
                 g, h, i};

  MatrixSymS<3,double> A(Adata, 6);

  UT_ASSERT_EQUALS(A(0,0), a);
  UT_ASSERT_EQUALS(A(0,1), b);
  UT_ASSERT_EQUALS(A(0,2), c);

  UT_ASSERT_EQUALS(A(1,0), d);
  UT_ASSERT_EQUALS(A(1,1), e);
  UT_ASSERT_EQUALS(A(1,2), f);

  UT_ASSERT_EQUALS(A(2,0), g);
  UT_ASSERT_EQUALS(A(2,1), h);
  UT_ASSERT_EQUALS(A(2,2), i);
}
UT_TEST_CASE_END( MatrixSymD_3x3_test )

//----------------------------------------------------------------------------//
UT_TEST_CASE( ATransposeA_ctor )
{
  const Int data[4] = {1,2,3,4};
  MatrixS<2,2,double> m1(data,4), m2(m1);

  MatrixSymS<2,double> sm1a;
  sm1a = Transpose(m1)*m1;

  UT_ASSERT_EQUALS(  10, sm1a(0,0) );
  UT_ASSERT_EQUALS(  14, sm1a(0,1) );
  UT_ASSERT_EQUALS(  14, sm1a(1,0) );
  UT_ASSERT_EQUALS(  20, sm1a(1,1) );

  sm1a += Transpose(m1)*m1;

  UT_ASSERT_EQUALS( 2*10, sm1a(0,0) );
  UT_ASSERT_EQUALS( 2*14, sm1a(0,1) );
  UT_ASSERT_EQUALS( 2*14, sm1a(1,0) );
  UT_ASSERT_EQUALS( 2*20, sm1a(1,1) );

  sm1a -= Transpose(m1)*m1;

  UT_ASSERT_EQUALS(  10, sm1a(0,0) );
  UT_ASSERT_EQUALS(  14, sm1a(0,1) );
  UT_ASSERT_EQUALS(  14, sm1a(1,0) );
  UT_ASSERT_EQUALS(  20, sm1a(1,1) );

  MatrixSymS<2,double> sm1b;
  sm1b = m1*Transpose(m1);

  UT_ASSERT_EQUALS(   5, sm1b(0,0) );
  UT_ASSERT_EQUALS(  11, sm1b(0,1) );
  UT_ASSERT_EQUALS(  11, sm1b(1,0) );
  UT_ASSERT_EQUALS(  25, sm1b(1,1) );

  sm1b += m1*Transpose(m1);

  UT_ASSERT_EQUALS(  2*5, sm1b(0,0) );
  UT_ASSERT_EQUALS( 2*11, sm1b(0,1) );
  UT_ASSERT_EQUALS( 2*11, sm1b(1,0) );
  UT_ASSERT_EQUALS( 2*25, sm1b(1,1) );

  sm1b -= m1*Transpose(m1);

  UT_ASSERT_EQUALS(   5, sm1b(0,0) );
  UT_ASSERT_EQUALS(  11, sm1b(0,1) );
  UT_ASSERT_EQUALS(  11, sm1b(1,0) );
  UT_ASSERT_EQUALS(  25, sm1b(1,1) );


  MatrixSymS<2,double> sm2a;
  sm2a = Transpose(m1)*m1;

  UT_ASSERT_EQUALS(  10, sm2a(0,0) );
  UT_ASSERT_EQUALS(  14, sm2a(0,1) );
  UT_ASSERT_EQUALS(  14, sm2a(1,0) );
  UT_ASSERT_EQUALS(  20, sm2a(1,1) );

  MatrixSymS<2,double> sm2b;
  sm2b = m1*Transpose(m1);

  UT_ASSERT_EQUALS(   5, sm2b(0,0) );
  UT_ASSERT_EQUALS(  11, sm2b(0,1) );
  UT_ASSERT_EQUALS(  11, sm2b(1,0) );
  UT_ASSERT_EQUALS(  25, sm2b(1,1) );


  VectorS<2,double> v1;
  v1 = {2,3};

  MatrixSymS<2,double> sm3a;
  sm3a = Transpose(m1)*diag(v1)*m1;

  UT_ASSERT_EQUALS(  29, sm3a(0,0) );
  UT_ASSERT_EQUALS(  40, sm3a(0,1) );
  UT_ASSERT_EQUALS(  40, sm3a(1,0) );
  UT_ASSERT_EQUALS(  56, sm3a(1,1) );

  sm3a += Transpose(m1)*diag(v1)*m1;

  UT_ASSERT_EQUALS( 2*29, sm3a(0,0) );
  UT_ASSERT_EQUALS( 2*40, sm3a(0,1) );
  UT_ASSERT_EQUALS( 2*40, sm3a(1,0) );
  UT_ASSERT_EQUALS( 2*56, sm3a(1,1) );

  sm3a -= Transpose(m1)*diag(v1)*m1;

  UT_ASSERT_EQUALS(  29, sm3a(0,0) );
  UT_ASSERT_EQUALS(  40, sm3a(0,1) );
  UT_ASSERT_EQUALS(  40, sm3a(1,0) );
  UT_ASSERT_EQUALS(  56, sm3a(1,1) );


  MatrixSymS<2,double> sm3b;
  sm3b = m1*diag(v1)*Transpose(m1);

  UT_ASSERT_EQUALS(  14, sm3b(0,0) );
  UT_ASSERT_EQUALS(  30, sm3b(0,1) );
  UT_ASSERT_EQUALS(  30, sm3b(1,0) );
  UT_ASSERT_EQUALS(  66, sm3b(1,1) );

  sm3b += m1*diag(v1)*Transpose(m1);

  UT_ASSERT_EQUALS( 2*14, sm3b(0,0) );
  UT_ASSERT_EQUALS( 2*30, sm3b(0,1) );
  UT_ASSERT_EQUALS( 2*30, sm3b(1,0) );
  UT_ASSERT_EQUALS( 2*66, sm3b(1,1) );

  sm3b -= m1*diag(v1)*Transpose(m1);

  UT_ASSERT_EQUALS(  14, sm3b(0,0) );
  UT_ASSERT_EQUALS(  30, sm3b(0,1) );
  UT_ASSERT_EQUALS(  30, sm3b(1,0) );
  UT_ASSERT_EQUALS(  66, sm3b(1,1) );


  MatrixSymS<2,double> sm4;
  sm4 = {{2},{4,3}};

  MatrixSymS<2,double> sm5a;
  sm5a = Transpose(m1)*sm4*m1;
  UT_ASSERT_EQUALS(  53, sm5a(0,0) );
  UT_ASSERT_EQUALS(  80, sm5a(0,1) );
  UT_ASSERT_EQUALS(  80, sm5a(1,0) );
  UT_ASSERT_EQUALS( 120, sm5a(1,1) );

  sm5a += Transpose(m1)*sm4*m1;
  UT_ASSERT_EQUALS(  2*53, sm5a(0,0) );
  UT_ASSERT_EQUALS(  2*80, sm5a(0,1) );
  UT_ASSERT_EQUALS(  2*80, sm5a(1,0) );
  UT_ASSERT_EQUALS( 2*120, sm5a(1,1) );

  sm5a -= Transpose(m1)*sm4*m1;
  UT_ASSERT_EQUALS(  53, sm5a(0,0) );
  UT_ASSERT_EQUALS(  80, sm5a(0,1) );
  UT_ASSERT_EQUALS(  80, sm5a(1,0) );
  UT_ASSERT_EQUALS( 120, sm5a(1,1) );


  MatrixSymS<2,double> sm5b;
  sm5b = m1*sm4*Transpose(m1);
  UT_ASSERT_EQUALS(  30, sm5b(0,0) );
  UT_ASSERT_EQUALS(  70, sm5b(0,1) );
  UT_ASSERT_EQUALS(  70, sm5b(1,0) );
  UT_ASSERT_EQUALS( 162, sm5b(1,1) );

  sm5b += m1*sm4*Transpose(m1);
  UT_ASSERT_EQUALS(  2*30, sm5b(0,0) );
  UT_ASSERT_EQUALS(  2*70, sm5b(0,1) );
  UT_ASSERT_EQUALS(  2*70, sm5b(1,0) );
  UT_ASSERT_EQUALS( 2*162, sm5b(1,1) );

  sm5b -= m1*sm4*Transpose(m1);
  UT_ASSERT_EQUALS(  30, sm5b(0,0) );
  UT_ASSERT_EQUALS(  70, sm5b(0,1) );
  UT_ASSERT_EQUALS(  70, sm5b(1,0) );
  UT_ASSERT_EQUALS( 162, sm5b(1,1) );


  sm1a = {{2},
          {3, 4}};

  sm1b = {{5},
          {6, 7}};

  MatrixSymS<2,double> sm6;
  sm6 = sm1a*sm1b*sm1a;

  UT_ASSERT_EQUALS(  155, sm6(0,0) );
  UT_ASSERT_EQUALS(  216, sm6(0,1) );
  UT_ASSERT_EQUALS(  216, sm6(1,0) );
  UT_ASSERT_EQUALS(  301, sm6(1,1) );

  sm6 += sm1a*sm1b*sm1a;

  UT_ASSERT_EQUALS( 2*155, sm6(0,0) );
  UT_ASSERT_EQUALS( 2*216, sm6(0,1) );
  UT_ASSERT_EQUALS( 2*216, sm6(1,0) );
  UT_ASSERT_EQUALS( 2*301, sm6(1,1) );

  sm6 -= sm1a*sm1b*sm1a;

  UT_ASSERT_EQUALS(  155, sm6(0,0) );
  UT_ASSERT_EQUALS(  216, sm6(0,1) );
  UT_ASSERT_EQUALS(  216, sm6(1,0) );
  UT_ASSERT_EQUALS(  301, sm6(1,1) );

  MatrixSymS<2,double> sm7;
  sm7 = sm1a*sm1a;

  UT_ASSERT_EQUALS( 13, sm7(0,0) );
  UT_ASSERT_EQUALS( 18, sm7(0,1) );
  UT_ASSERT_EQUALS( 18, sm7(1,0) );
  UT_ASSERT_EQUALS( 25, sm7(1,1) );

  sm7 += sm1a*sm1a;

  UT_ASSERT_EQUALS( 2*13, sm7(0,0) );
  UT_ASSERT_EQUALS( 2*18, sm7(0,1) );
  UT_ASSERT_EQUALS( 2*18, sm7(1,0) );
  UT_ASSERT_EQUALS( 2*25, sm7(1,1) );

  sm7 -= sm1a*sm1a;

  UT_ASSERT_EQUALS( 13, sm7(0,0) );
  UT_ASSERT_EQUALS( 18, sm7(0,1) );
  UT_ASSERT_EQUALS( 18, sm7(1,0) );
  UT_ASSERT_EQUALS( 25, sm7(1,1) );


  //Make sure you can't construct a potentially non-symmetric matrix
  /*UT_ASSERT_THROW( sm1a  = Transpose(m1)*m2, AssertionException );
  UT_ASSERT_THROW( sm1a += Transpose(m1)*m2, AssertionException );
  UT_ASSERT_THROW( sm1a -= Transpose(m1)*m2, AssertionException );

  UT_ASSERT_THROW( sm1b  = Transpose(m1)*diag(v1)*m2, AssertionException );
  UT_ASSERT_THROW( sm1b += Transpose(m1)*diag(v1)*m2, AssertionException );
  UT_ASSERT_THROW( sm1b -= Transpose(m1)*diag(v1)*m2, AssertionException );

  UT_ASSERT_THROW( sm1b  = Transpose(m1)*sm4*m2, AssertionException );
  UT_ASSERT_THROW( sm1b += Transpose(m1)*sm4*m2, AssertionException );
  UT_ASSERT_THROW( sm1b -= Transpose(m1)*sm4*m2, AssertionException );

  UT_ASSERT_THROW( sm1b  = sm1a*sm1b*sm2a, AssertionException );
  UT_ASSERT_THROW( sm1b += sm1a*sm1b*sm2a, AssertionException );
  UT_ASSERT_THROW( sm1b -= sm1a*sm1b*sm2a, AssertionException );

  UT_ASSERT_THROW( sm1b  = sm1a*sm1b, AssertionException );
  UT_ASSERT_THROW( sm1b += sm1a*sm1b, AssertionException );
  UT_ASSERT_THROW( sm1b -= sm1a*sm1b, AssertionException );*/
}
UT_TEST_CASE_END( ATransposeA_ctor )


#if 0



//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( ATransposeA_rectangle_ctor )
{

  const Int data[6] = {2,3,4,
                       5,6,7};
  MatrixS<2,3,Int> ma(data,6);

  MatrixSymS<2,Int> sm1a;
  sm1a = ma*Transpose(ma);
  UT_ASSERT_EQUALS(  29, sm1a(0,0) );
  UT_ASSERT_EQUALS(  56, sm1a(0,1) );
  UT_ASSERT_EQUALS(  56, sm1a(1,0) );
  UT_ASSERT_EQUALS( 110, sm1a(1,1) );

  sm1a += ma*Transpose(ma);

  UT_ASSERT_EQUALS(  2*29, sm1a(0,0) );
  UT_ASSERT_EQUALS(  2*56, sm1a(0,1) );
  UT_ASSERT_EQUALS(  2*56, sm1a(1,0) );
  UT_ASSERT_EQUALS( 2*110, sm1a(1,1) );

  sm1a -= ma*Transpose(ma);

  UT_ASSERT_EQUALS(  29, sm1a(0,0) );
  UT_ASSERT_EQUALS(  56, sm1a(0,1) );
  UT_ASSERT_EQUALS(  56, sm1a(1,0) );
  UT_ASSERT_EQUALS( 110, sm1a(1,1) );

  MatrixSymS<3,Int> sm1b;
  sm1b = Transpose(ma)*ma;

  UT_ASSERT_EQUALS( 29, sm1b(0,0) );
  UT_ASSERT_EQUALS( 36, sm1b(0,1) );
  UT_ASSERT_EQUALS( 43, sm1b(0,2) );
  UT_ASSERT_EQUALS( 36, sm1b(1,0) );
  UT_ASSERT_EQUALS( 45, sm1b(1,1) );
  UT_ASSERT_EQUALS( 54, sm1b(1,2) );
  UT_ASSERT_EQUALS( 43, sm1b(2,0) );
  UT_ASSERT_EQUALS( 54, sm1b(2,1) );
  UT_ASSERT_EQUALS( 65, sm1b(2,2) );

  sm1b += Transpose(ma)*ma;

  UT_ASSERT_EQUALS( 2*29, sm1b(0,0) );
  UT_ASSERT_EQUALS( 2*36, sm1b(0,1) );
  UT_ASSERT_EQUALS( 2*43, sm1b(0,2) );
  UT_ASSERT_EQUALS( 2*36, sm1b(1,0) );
  UT_ASSERT_EQUALS( 2*45, sm1b(1,1) );
  UT_ASSERT_EQUALS( 2*54, sm1b(1,2) );
  UT_ASSERT_EQUALS( 2*43, sm1b(2,0) );
  UT_ASSERT_EQUALS( 2*54, sm1b(2,1) );
  UT_ASSERT_EQUALS( 2*65, sm1b(2,2) );

  sm1b -= Transpose(ma)*ma;

  UT_ASSERT_EQUALS( 29, sm1b(0,0) );
  UT_ASSERT_EQUALS( 36, sm1b(0,1) );
  UT_ASSERT_EQUALS( 43, sm1b(0,2) );
  UT_ASSERT_EQUALS( 36, sm1b(1,0) );
  UT_ASSERT_EQUALS( 45, sm1b(1,1) );
  UT_ASSERT_EQUALS( 54, sm1b(1,2) );
  UT_ASSERT_EQUALS( 43, sm1b(2,0) );
  UT_ASSERT_EQUALS( 54, sm1b(2,1) );
  UT_ASSERT_EQUALS( 65, sm1b(2,2) );

  const Int data2[6] = {7,6,
                        5,4,
                        3,2};
  MatrixS<3,2,Int> mb(data2,6);
  MatrixSymS<2,Int> sm2a;

  sm2a = Transpose(mb)*mb;

  UT_ASSERT_EQUALS( 83, sm2a(0,0) );
  UT_ASSERT_EQUALS( 68, sm2a(0,1) );
  UT_ASSERT_EQUALS( 68, sm2a(1,0) );
  UT_ASSERT_EQUALS( 56, sm2a(1,1) );

  sm2a += Transpose(mb)*mb;

  UT_ASSERT_EQUALS( 2*83, sm2a(0,0) );
  UT_ASSERT_EQUALS( 2*68, sm2a(0,1) );
  UT_ASSERT_EQUALS( 2*68, sm2a(1,0) );
  UT_ASSERT_EQUALS( 2*56, sm2a(1,1) );

  sm2a -= Transpose(mb)*mb;

  UT_ASSERT_EQUALS( 83, sm2a(0,0) );
  UT_ASSERT_EQUALS( 68, sm2a(0,1) );
  UT_ASSERT_EQUALS( 68, sm2a(1,0) );
  UT_ASSERT_EQUALS( 56, sm2a(1,1) );

  MatrixSymS<3,Int> sm2b;

  sm2b = mb*Transpose(mb);

  UT_ASSERT_EQUALS( 85, sm2b(0,0) );
  UT_ASSERT_EQUALS( 59, sm2b(0,1) );
  UT_ASSERT_EQUALS( 33, sm2b(0,2) );
  UT_ASSERT_EQUALS( 59, sm2b(1,0) );
  UT_ASSERT_EQUALS( 41, sm2b(1,1) );
  UT_ASSERT_EQUALS( 23, sm2b(1,2) );
  UT_ASSERT_EQUALS( 33, sm2b(2,0) );
  UT_ASSERT_EQUALS( 23, sm2b(2,1) );
  UT_ASSERT_EQUALS( 13, sm2b(2,2) );

  sm2b += mb*Transpose(mb);

  UT_ASSERT_EQUALS( 2*85, sm2b(0,0) );
  UT_ASSERT_EQUALS( 2*59, sm2b(0,1) );
  UT_ASSERT_EQUALS( 2*33, sm2b(0,2) );
  UT_ASSERT_EQUALS( 2*59, sm2b(1,0) );
  UT_ASSERT_EQUALS( 2*41, sm2b(1,1) );
  UT_ASSERT_EQUALS( 2*23, sm2b(1,2) );
  UT_ASSERT_EQUALS( 2*33, sm2b(2,0) );
  UT_ASSERT_EQUALS( 2*23, sm2b(2,1) );
  UT_ASSERT_EQUALS( 2*13, sm2b(2,2) );

  sm2b -= mb*Transpose(mb);

  UT_ASSERT_EQUALS( 85, sm2b(0,0) );
  UT_ASSERT_EQUALS( 59, sm2b(0,1) );
  UT_ASSERT_EQUALS( 33, sm2b(0,2) );
  UT_ASSERT_EQUALS( 59, sm2b(1,0) );
  UT_ASSERT_EQUALS( 41, sm2b(1,1) );
  UT_ASSERT_EQUALS( 23, sm2b(1,2) );
  UT_ASSERT_EQUALS( 33, sm2b(2,0) );
  UT_ASSERT_EQUALS( 23, sm2b(2,1) );
  UT_ASSERT_EQUALS( 13, sm2b(2,2) );

  VectorS2 v1 = {1,1};
  VectorS<3,Int> v2 = {1,1,1};
  MatrixSymS2 sm3a;

  sm3a = ma*diag(v2)*Transpose(ma);

  UT_ASSERT_EQUALS(  29, sm3a(0,0) );
  UT_ASSERT_EQUALS(  56, sm3a(0,1) );
  UT_ASSERT_EQUALS(  56, sm3a(1,0) );
  UT_ASSERT_EQUALS( 110, sm3a(1,1) );

  sm3a += ma*diag(v2)*Transpose(ma);

  UT_ASSERT_EQUALS(  2*29, sm3a(0,0) );
  UT_ASSERT_EQUALS(  2*56, sm3a(0,1) );
  UT_ASSERT_EQUALS(  2*56, sm3a(1,0) );
  UT_ASSERT_EQUALS( 2*110, sm3a(1,1) );

  sm3a -= ma*diag(v2)*Transpose(ma);

  UT_ASSERT_EQUALS(  29, sm3a(0,0) );
  UT_ASSERT_EQUALS(  56, sm3a(0,1) );
  UT_ASSERT_EQUALS(  56, sm3a(1,0) );
  UT_ASSERT_EQUALS( 110, sm3a(1,1) );

  MatrixSymS<3,Int> sm3b;
  sm3b = Transpose(ma)*diag(v1)*ma;

  UT_ASSERT_EQUALS( 29, sm3b(0,0) );
  UT_ASSERT_EQUALS( 36, sm3b(0,1) );
  UT_ASSERT_EQUALS( 43, sm3b(0,2) );
  UT_ASSERT_EQUALS( 36, sm3b(1,0) );
  UT_ASSERT_EQUALS( 45, sm3b(1,1) );
  UT_ASSERT_EQUALS( 54, sm3b(1,2) );
  UT_ASSERT_EQUALS( 43, sm3b(2,0) );
  UT_ASSERT_EQUALS( 54, sm3b(2,1) );
  UT_ASSERT_EQUALS( 65, sm3b(2,2) );

  sm3b += Transpose(ma)*diag(v1)*ma;

  UT_ASSERT_EQUALS( 2*29, sm3b(0,0) );
  UT_ASSERT_EQUALS( 2*36, sm3b(0,1) );
  UT_ASSERT_EQUALS( 2*43, sm3b(0,2) );
  UT_ASSERT_EQUALS( 2*36, sm3b(1,0) );
  UT_ASSERT_EQUALS( 2*45, sm3b(1,1) );
  UT_ASSERT_EQUALS( 2*54, sm3b(1,2) );
  UT_ASSERT_EQUALS( 2*43, sm3b(2,0) );
  UT_ASSERT_EQUALS( 2*54, sm3b(2,1) );
  UT_ASSERT_EQUALS( 2*65, sm3b(2,2) );

  sm3b -= Transpose(ma)*diag(v1)*ma;

  UT_ASSERT_EQUALS( 29, sm3b(0,0) );
  UT_ASSERT_EQUALS( 36, sm3b(0,1) );
  UT_ASSERT_EQUALS( 43, sm3b(0,2) );
  UT_ASSERT_EQUALS( 36, sm3b(1,0) );
  UT_ASSERT_EQUALS( 45, sm3b(1,1) );
  UT_ASSERT_EQUALS( 54, sm3b(1,2) );
  UT_ASSERT_EQUALS( 43, sm3b(2,0) );
  UT_ASSERT_EQUALS( 54, sm3b(2,1) );
  UT_ASSERT_EQUALS( 65, sm3b(2,2) );

  MatrixSymS<2,Int> sm4a;

  sm4a = Transpose(mb)*diag(v2)*mb;

  UT_ASSERT_EQUALS( 83, sm2a(0,0) );
  UT_ASSERT_EQUALS( 68, sm2a(0,1) );
  UT_ASSERT_EQUALS( 68, sm2a(1,0) );
  UT_ASSERT_EQUALS( 56, sm2a(1,1) );

  sm4a += Transpose(mb)*diag(v2)*mb;

  UT_ASSERT_EQUALS( 2*83, sm4a(0,0) );
  UT_ASSERT_EQUALS( 2*68, sm4a(0,1) );
  UT_ASSERT_EQUALS( 2*68, sm4a(1,0) );
  UT_ASSERT_EQUALS( 2*56, sm4a(1,1) );

  sm4a -= Transpose(mb)*diag(v2)*mb;

  UT_ASSERT_EQUALS( 83, sm2a(0,0) );
  UT_ASSERT_EQUALS( 68, sm2a(0,1) );
  UT_ASSERT_EQUALS( 68, sm2a(1,0) );
  UT_ASSERT_EQUALS( 56, sm2a(1,1) );

  MatrixSymS<3,Int> sm4b;

  sm4b = mb*diag(v1)*Transpose(mb);

  UT_ASSERT_EQUALS( 85, sm4b(0,0) );
  UT_ASSERT_EQUALS( 59, sm4b(0,1) );
  UT_ASSERT_EQUALS( 33, sm4b(0,2) );
  UT_ASSERT_EQUALS( 59, sm4b(1,0) );
  UT_ASSERT_EQUALS( 41, sm4b(1,1) );
  UT_ASSERT_EQUALS( 23, sm4b(1,2) );
  UT_ASSERT_EQUALS( 33, sm4b(2,0) );
  UT_ASSERT_EQUALS( 23, sm4b(2,1) );
  UT_ASSERT_EQUALS( 13, sm4b(2,2) );

  sm4b += mb*diag(v1)*Transpose(mb);

  UT_ASSERT_EQUALS( 2*85, sm4b(0,0) );
  UT_ASSERT_EQUALS( 2*59, sm4b(0,1) );
  UT_ASSERT_EQUALS( 2*33, sm4b(0,2) );
  UT_ASSERT_EQUALS( 2*59, sm4b(1,0) );
  UT_ASSERT_EQUALS( 2*41, sm4b(1,1) );
  UT_ASSERT_EQUALS( 2*23, sm4b(1,2) );
  UT_ASSERT_EQUALS( 2*33, sm4b(2,0) );
  UT_ASSERT_EQUALS( 2*23, sm4b(2,1) );
  UT_ASSERT_EQUALS( 2*13, sm4b(2,2) );

  sm4b -= mb*diag(v1)*Transpose(mb);

  UT_ASSERT_EQUALS( 85, sm4b(0,0) );
  UT_ASSERT_EQUALS( 59, sm4b(0,1) );
  UT_ASSERT_EQUALS( 33, sm4b(0,2) );
  UT_ASSERT_EQUALS( 59, sm4b(1,0) );
  UT_ASSERT_EQUALS( 41, sm4b(1,1) );
  UT_ASSERT_EQUALS( 23, sm4b(1,2) );
  UT_ASSERT_EQUALS( 33, sm4b(2,0) );
  UT_ASSERT_EQUALS( 23, sm4b(2,1) );
  UT_ASSERT_EQUALS( 13, sm4b(2,2) );

  MatrixSymS2 sma = {{2},{4,3}};

  MatrixSymS<3,Int> sm5a;

  sm5a = Transpose(ma)*sma*ma;

  UT_ASSERT_EQUALS( 163, sm5a(0,0) );
  UT_ASSERT_EQUALS( 210, sm5a(0,1) );
  UT_ASSERT_EQUALS( 257, sm5a(0,2) );
  UT_ASSERT_EQUALS( 210, sm5a(1,0) );
  UT_ASSERT_EQUALS( 270, sm5a(1,1) );
  UT_ASSERT_EQUALS( 330, sm5a(1,2) );
  UT_ASSERT_EQUALS( 257, sm5a(2,0) );
  UT_ASSERT_EQUALS( 330, sm5a(2,1) );
  UT_ASSERT_EQUALS( 403, sm5a(2,2) );

  sm5a += Transpose(ma)*sma*ma;

  UT_ASSERT_EQUALS( 2*163, sm5a(0,0) );
  UT_ASSERT_EQUALS( 2*210, sm5a(0,1) );
  UT_ASSERT_EQUALS( 2*257, sm5a(0,2) );
  UT_ASSERT_EQUALS( 2*210, sm5a(1,0) );
  UT_ASSERT_EQUALS( 2*270, sm5a(1,1) );
  UT_ASSERT_EQUALS( 2*330, sm5a(1,2) );
  UT_ASSERT_EQUALS( 2*257, sm5a(2,0) );
  UT_ASSERT_EQUALS( 2*330, sm5a(2,1) );
  UT_ASSERT_EQUALS( 2*403, sm5a(2,2) );

  sm5a -= Transpose(ma)*sma*ma;

  UT_ASSERT_EQUALS( 163, sm5a(0,0) );
  UT_ASSERT_EQUALS( 210, sm5a(0,1) );
  UT_ASSERT_EQUALS( 257, sm5a(0,2) );
  UT_ASSERT_EQUALS( 210, sm5a(1,0) );
  UT_ASSERT_EQUALS( 270, sm5a(1,1) );
  UT_ASSERT_EQUALS( 330, sm5a(1,2) );
  UT_ASSERT_EQUALS( 257, sm5a(2,0) );
  UT_ASSERT_EQUALS( 330, sm5a(2,1) );
  UT_ASSERT_EQUALS( 403, sm5a(2,2) );

  MatrixSymS<3,Int> smb = {{2},
                           {4,3},
                           {4,5,8}};

  MatrixSymS<2,Int> sm5b;

  sm5b = ma*smb*Transpose(ma);

  UT_ASSERT_EQUALS(  395, sm5b(0,0) );
  UT_ASSERT_EQUALS(  767, sm5b(0,1) );
  UT_ASSERT_EQUALS(  767, sm5b(1,0) );
  UT_ASSERT_EQUALS( 1490, sm5b(1,1) );

  sm5b += ma*smb*Transpose(ma);

  UT_ASSERT_EQUALS(  2*395, sm5b(0,0) );
  UT_ASSERT_EQUALS(  2*767, sm5b(0,1) );
  UT_ASSERT_EQUALS(  2*767, sm5b(1,0) );
  UT_ASSERT_EQUALS( 2*1490, sm5b(1,1) );

  sm5b -= ma*smb*Transpose(ma);

  UT_ASSERT_EQUALS(  395, sm5b(0,0) );
  UT_ASSERT_EQUALS(  767, sm5b(0,1) );
  UT_ASSERT_EQUALS(  767, sm5b(1,0) );
  UT_ASSERT_EQUALS( 1490, sm5b(1,1) );

  MatrixSymS<3,Int> sm6a;

  sm6a = mb*sma*Transpose(mb);

  UT_ASSERT_EQUALS( 542, sm6a(0,0) );
  UT_ASSERT_EQUALS( 374, sm6a(0,1) );
  UT_ASSERT_EQUALS( 206, sm6a(0,2) );
  UT_ASSERT_EQUALS( 374, sm6a(1,0) );
  UT_ASSERT_EQUALS( 258, sm6a(1,1) );
  UT_ASSERT_EQUALS( 142, sm6a(1,2) );
  UT_ASSERT_EQUALS( 206, sm6a(2,0) );
  UT_ASSERT_EQUALS( 142, sm6a(2,1) );
  UT_ASSERT_EQUALS(  78, sm6a(2,2) );

  sm6a += mb*sma*Transpose(mb);

  UT_ASSERT_EQUALS( 2*542, sm6a(0,0) );
  UT_ASSERT_EQUALS( 2*374, sm6a(0,1) );
  UT_ASSERT_EQUALS( 2*206, sm6a(0,2) );
  UT_ASSERT_EQUALS( 2*374, sm6a(1,0) );
  UT_ASSERT_EQUALS( 2*258, sm6a(1,1) );
  UT_ASSERT_EQUALS( 2*142, sm6a(1,2) );
  UT_ASSERT_EQUALS( 2*206, sm6a(2,0) );
  UT_ASSERT_EQUALS( 2*142, sm6a(2,1) );
  UT_ASSERT_EQUALS(  2*78, sm6a(2,2) );

  sm6a -= mb*sma*Transpose(mb);

  UT_ASSERT_EQUALS( 542, sm6a(0,0) );
  UT_ASSERT_EQUALS( 374, sm6a(0,1) );
  UT_ASSERT_EQUALS( 206, sm6a(0,2) );
  UT_ASSERT_EQUALS( 374, sm6a(1,0) );
  UT_ASSERT_EQUALS( 258, sm6a(1,1) );
  UT_ASSERT_EQUALS( 142, sm6a(1,2) );
  UT_ASSERT_EQUALS( 206, sm6a(2,0) );
  UT_ASSERT_EQUALS( 142, sm6a(2,1) );
  UT_ASSERT_EQUALS(  78, sm6a(2,2) );

  MatrixSymS<2,Int> sm6b;

  sm6b = Transpose(mb)*smb*mb;

  UT_ASSERT_EQUALS( 843, sm6b(0,0) );
  UT_ASSERT_EQUALS( 662, sm6b(0,1) );
  UT_ASSERT_EQUALS( 662, sm6b(1,0) );
  UT_ASSERT_EQUALS( 520, sm6b(1,1) );

  sm6b += Transpose(mb)*smb*mb;

  UT_ASSERT_EQUALS( 2*843, sm6b(0,0) );
  UT_ASSERT_EQUALS( 2*662, sm6b(0,1) );
  UT_ASSERT_EQUALS( 2*662, sm6b(1,0) );
  UT_ASSERT_EQUALS( 2*520, sm6b(1,1) );

  sm6b -= Transpose(mb)*smb*mb;

  UT_ASSERT_EQUALS( 843, sm6b(0,0) );
  UT_ASSERT_EQUALS( 662, sm6b(0,1) );
  UT_ASSERT_EQUALS( 662, sm6b(1,0) );
  UT_ASSERT_EQUALS( 520, sm6b(1,1) );

  MatrixS<3,2,Int> mc(data,6);
  MatrixS<2,3,Int> md(data,6);

  //Make sure you can't construct a potentially non-symmetric matrix
  UT_ASSERT_THROW( sm1a  = Transpose(mb)*mc;, AssertionException );
  UT_ASSERT_THROW( sm1a += Transpose(mb)*mc;, AssertionException );
  UT_ASSERT_THROW( sm1a -= Transpose(mb)*mc;, AssertionException );

  UT_ASSERT_THROW( sm1b  = Transpose(md)*diag(v1)*ma;, AssertionException );
  UT_ASSERT_THROW( sm1b += Transpose(md)*diag(v1)*ma;, AssertionException );
  UT_ASSERT_THROW( sm1b -= Transpose(md)*diag(v1)*ma;, AssertionException );

  UT_ASSERT_THROW( sm1b  = Transpose(md)*sm6b*ma;, AssertionException );
  UT_ASSERT_THROW( sm1b += Transpose(md)*sm6b*ma;, AssertionException );
  UT_ASSERT_THROW( sm1b -= Transpose(md)*sm6b*ma;, AssertionException );
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

  UT_ASSERT_EQUALS( xTMx , 145 );
}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( matrix_ops2 )
{
  MatrixSymS2 m1 = {{1},
                    {2,4}};
  MatrixSymS2 m2(m1), m3, m4, m5;

  // size
  UT_ASSERT( m1.M == 2 );
  UT_ASSERT( m1.N == 2 );
  UT_ASSERT( m2.M == 2 );
  UT_ASSERT( m2.N == 2 );

  // ctors
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );

  // assignment
  m3 = m1;
  m4 = 5;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 5,5,5,5 ) );

  m2 = m3 = 3;
  UT_ASSERT( chkMatrixSymS22( m2, 3,3,3,3 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 3,3,3,3 ) );

  // unary
  m2 = +m1;
  m3 = -m1;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );

  // binary accumulation
  m4 = m1;
  m4 *= 5;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 5,10,10,20 ) );

  m2 = 5;
  m3 = 5;
  m2 += m1;
  m3 -= m1;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 6,7,7,9 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 4,3,3,1 ) );

  // binary operators
//  m2 = m1 + 3;
//  m3 = m1 - 3;
  m4 = m1 * 3;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
//  UT_ASSERT( chkMatrixSymS22( m2, 4,5,6,7 ) );
//  UT_ASSERT( chkMatrixSymS22( m3, -2,-1,0,1 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 3,6,6,12 ) );

//  m2 = 3 + m1;
//  m3 = 3 - m1;
  m4 = 3 * m1;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
//  UT_ASSERT( chkMatrixSymS22( m2, 4,5,6,7 ) );
//  UT_ASSERT( chkMatrixSymS22( m3, 2,1,0,-1 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 3,6,6,12 ) );

  m2 = 3;
  m3 = m1 + m2;
  m4 = m1 - m2;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 3,3,3,3 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 4,5,5,7 ) );
  UT_ASSERT( chkMatrixSymS22( m4, -2,-1,-1,1 ) );

  // arithmetic combinations

  m2 = m1;
  m3 = m1 + m2;
  m4 = m1 + m2 + m3;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 2,4,4,8 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 4,8,8,16 ) );

  m2 += m1;
  m3 += m1 + m2;
  m4 += m1 + m2 + m3;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 2,4,4,8 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 5,10,10,20 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 12,24,24,48 ) );

  m3 = m1 - m2;
  m4 = m1 - m2 - m3;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 2,4,4,8 ) );
  UT_ASSERT( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 0,0,0,0 ) );

  m2 -= m1;
  m3 -= m1 - m2;
  m4 -= m1 - m2 - m3;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );
  UT_ASSERT( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );

  m3 = m1 - m2;
  m4 = m1 + m2 - m3;
  m5 = m1 - m2 + m3;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );
  UT_ASSERT( chkMatrixSymS22( m5, 0,0,0,0 ) );

  m5 = (m1 + m2) + (m3 + m4);
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = (m1 + m2) + (m3 - m4);
  UT_ASSERT( chkMatrixSymS22( m5, 0,0,0,0 ) );
  m5 = (m1 + m2) - (m3 + m4);
  UT_ASSERT( chkMatrixSymS22( m5, 0,0,0,0 ) );
  m5 = (m1 + m2) - (m3 - m4);
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = (m1 - m2) + (m3 + m4);
  UT_ASSERT( chkMatrixSymS22( m5, 2,4,4,8 ) );
  m5 = (m1 - m2) + (m3 - m4);
  UT_ASSERT( chkMatrixSymS22( m5, -2,-4,-4,-8 ) );
  m5 = (m1 - m2) - (m3 + m4);
  UT_ASSERT( chkMatrixSymS22( m5, -2,-4,-4,-8 ) );
  m5 = (m1 - m2) - (m3 - m4);
  UT_ASSERT( chkMatrixSymS22( m5, 2,4,4,8 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 += (m1 + m2) + (m3 + m4);
  m5 += (m1 + m2) + (m3 - m4);
  m5 += (m1 + m2) - (m3 + m4);
  m5 += (m1 + m2) - (m3 - m4);
  m5 += (m1 - m2) + (m3 + m4);
  m5 += (m1 - m2) + (m3 - m4);
  m5 += (m1 - m2) - (m3 + m4);
  m5 += (m1 - m2) - (m3 - m4);
  UT_ASSERT( chkMatrixSymS22( m5, 10,20,20,40 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 -= (m1 + m2) + (m3 + m4);
  m5 -= (m1 + m2) + (m3 - m4);
  m5 -= (m1 + m2) - (m3 + m4);
  m5 -= (m1 + m2) - (m3 - m4);
  m5 -= (m1 - m2) + (m3 + m4);
  m5 -= (m1 - m2) + (m3 - m4);
  m5 -= (m1 - m2) - (m3 + m4);
  m5 -= (m1 - m2) - (m3 - m4);
  UT_ASSERT( chkMatrixSymS22( m5, 2,4,4,8 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 0,0,0,0 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m2 = 1*m1;
  m3 = m2*2;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 2,4,4,8 ) );

  m2 += 1*m1;
  m3 += m2*2;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 2,4,4,8 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 6,12,12,24 ) );

  m2 -= 1*m1;
  m3 -= m2*2;
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 4,8,8,16 ) );



  m3 = (2*m2)*2;
  UT_ASSERT( chkMatrixSymS22( m3, 4,8,8,16 ) );

  m3 += (2*m2)*2;
  UT_ASSERT( chkMatrixSymS22( m3, 2*4,2*8,2*8,2*16 ) );

  m3 = 2*(m2*2);
  UT_ASSERT( chkMatrixSymS22( m3, 4,8,8,16 ) );

  m3 += 2*(m2*2);
  UT_ASSERT( chkMatrixSymS22( m3, 2*4,2*8,2*8,2*16 ) );

  m3 = (4*m2)/2;
  UT_ASSERT( chkMatrixSymS22( m3, 2,4,4,8 ) );

  m3 += (4*m2)/2;
  UT_ASSERT( chkMatrixSymS22( m3, 2*2,2*4,2*4,2*8 ) );



  m5 = 1*(m1 + m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 14,28,28,56 ) );
  m5 = 1*(m1 + m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 6,12,12,24 ) );
  m5 = 1*(m1 + m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -10,-20,-20,-40 ) );
  m5 = 1*(m1 + m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -2,-4,-4,-8 ) );
  m5 = 1*(m1 - m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 12,24,24,48 ) );
  m5 = 1*(m1 - m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = 1*(m1 - m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -12,-24,-24,-48 ) );
  m5 = 1*(m1 - m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -4,-8,-8,-16 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 4,8,8,16 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 += 1*(m1 + m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 10,20,20,40 ) );
  m5 += 1*(m1 + m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 16,32,32,64 ) );
  m5 += 1*(m1 + m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 6,12,12,24 ) );
  m5 += 1*(m1 + m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 += 1*(m1 - m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 16,32,32,64 ) );
  m5 += 1*(m1 - m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 20,40,40,80 ) );
  m5 += 1*(m1 - m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 8,16,16,32 ) );
  m5 += 1*(m1 - m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 4,8,8,16 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 -= 1*(m1 + m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -10,-20,-20,-40 ) );
  m5 -= 1*(m1 + m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -16,-32,-32,-64 ) );
  m5 -= 1*(m1 + m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -6,-12,-12,-24 ) );
  m5 -= 1*(m1 + m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -4,-8,-8,-16 ) );
  m5 -= 1*(m1 - m2) + (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -16,-32,-32,-64 ) );
  m5 -= 1*(m1 - m2) + (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -20,-40,-40,-80 ) );
  m5 -= 1*(m1 - m2) - (m3 + m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -8,-16,-16,-32 ) );
  m5 -= 1*(m1 - m2) - (m3 - m4)*2;
  UT_ASSERT( chkMatrixSymS22( m5, -4,-8,-8,-16 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, 4,8,8,16 ) );
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );

  m5 = 1*(m1 + m2)*2;
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = 1*2*(m1 + m2);
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  m5 = (m1 + m2)*1*2;
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );

  //Need to divide by 1 because we are working with integers
  m5 = (m1 + m2)/1;
  UT_ASSERT( chkMatrixSymS22( m5, 2,4,4,8 ) );

  m5 = +( 2*(m1 + m2) );
  UT_ASSERT( chkMatrixSymS22( m5, 4,8,8,16 ) );

  m2 = +m1;
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  m3 = -m2;
  UT_ASSERT( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );
  m4 = +(m1 + m2);
  UT_ASSERT( chkMatrixSymS22( m4, 2,4,4,8 ) );
  m4 = +(m1 - m2);
  UT_ASSERT( chkMatrixSymS22( m4, 0,0,0,0 ) );
  m4 = -(m1 + m2);
  UT_ASSERT( chkMatrixSymS22( m4, -2,-4,-4,-8 ) );
  m4 = -(m1 - m2);
  UT_ASSERT( chkMatrixSymS22( m4, 0,0,0,0 ) );
  m4 = +(m1 + m2) + m3;
  UT_ASSERT( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -(m1 + m2) + m3;
  UT_ASSERT( chkMatrixSymS22( m4, -3,-6,-6,-12 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 1,2,2,4 ) );
  UT_ASSERT( chkMatrixSymS22( m3, -1,-2,-2,-4 ) );

  m4 = +1*m1;
  UT_ASSERT( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -1*m1;
  UT_ASSERT( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  m4 = +m1*1;
  UT_ASSERT( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -m1*1;
  UT_ASSERT( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  m4 = +(1*m1);
  UT_ASSERT( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -(1*m1);
  UT_ASSERT( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  m4 = +(m1*1);
  UT_ASSERT( chkMatrixSymS22( m4, 1,2,2,4 ) );
  m4 = -(m1*1);
  UT_ASSERT( chkMatrixSymS22( m4, -1,-2,-2,-4 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 1,2,2,4 ) );
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

  UT_ASSERT( chkVectorS2( a1, 1,2 ) );
  UT_ASSERT( chkVectorS2( a2, 1,2 ) );

  UT_ASSERT( chkMatrixSymS22( m1, 3,4,4,6 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 3,4,4,6 ) );

  a2 = m1*a1;
  UT_ASSERT( chkVectorS2( a2, 11,16 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 3,4,4,6 ) );

  a2 = +(m1*a1);
  UT_ASSERT( chkVectorS2( a2, 11,16 ) );

  a2 += m1*a1;
  UT_ASSERT( chkVectorS2( a2, 2*11,2*16 ) );

  a2 -= m1*a1;
  UT_ASSERT( chkVectorS2( a2, 11,16 ) );

  a2 += m1*a1 + m1*a1;
  UT_ASSERT( chkVectorS2( a2, 3*11,3*16 ) );

  a2 += m1*a1 - m1*a1;
  UT_ASSERT( chkVectorS2( a2, 3*11,3*16 ) );



  a2 = (2*m1)*a1;
  UT_ASSERT( chkVectorS2( a2, 2*11,2*16 ) );

  a2 = +((2*m1)*a1);
  UT_ASSERT( chkVectorS2( a2, 2*11,2*16 ) );

  a2 += (2*m1)*a1;
  UT_ASSERT( chkVectorS2( a2, 4*11,4*16 ) );

  a2 = m1*(a1*2);
  UT_ASSERT( chkVectorS2( a2, 2*11,2*16 ) );

  a2 = +(m1*(a1*2));
  UT_ASSERT( chkVectorS2( a2, 2*11,2*16 ) );

  a2 += m1*(a1*2);
  UT_ASSERT( chkVectorS2( a2, 4*11,4*16 ) );



  a2 = (m1 + m2)*a1;
  UT_ASSERT( chkVectorS2( a2, 22,32 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 3,4,4,6 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 3,4,4,6 ) );

  a2 = +((m1 + m2)*a1);
  UT_ASSERT( chkVectorS2( a2, 22,32 ) );

  a2 += (m1 + m2)*a1;
  UT_ASSERT( chkVectorS2( a2, 2*22,2*32 ) );

  a2 -= (m1 + m2)*a1;
  UT_ASSERT( chkVectorS2( a2, 22,32 ) );

//  a2 = a1*(m1 + m2);
//  UT_ASSERT( chkVectorS2( a2, 26,32 ) );
//  UT_ASSERT( chkMatrixSymS22( m1, 3,4,5,6 ) );
//  UT_ASSERT( chkMatrixSymS22( m2, 3,4,5,6 ) );

  a2 = (m1 - m2)*a1;
  UT_ASSERT( chkVectorS2( a2, 0,0 ) );
  UT_ASSERT( chkMatrixSymS22( m1, 3,4,4,6 ) );
  UT_ASSERT( chkMatrixSymS22( m2, 3,4,4,6 ) );

//  a2 = a1*(m1 - m2);
//  UT_ASSERT( chkVectorS2( a2, 0,0 ) );
//  UT_ASSERT( chkMatrixSymS22( m1, 3,4,5,6 ) );
//  UT_ASSERT( chkMatrixSymS22( m2, 3,4,5,6 ) );

  a2 = a1;

  a3 = m1*(a1 + a2);
  UT_ASSERT( chkVectorS2( a3, 2*11,2*16 ) );

  a3 = +(m1*(a1 + a2));
  UT_ASSERT( chkVectorS2( a3, 2*11,2*16 ) );

  a3 += m1*(a1 + a2);
  UT_ASSERT( chkVectorS2( a3, 4*11,4*16 ) );

  a3 -= m1*(a1 + a2);
  UT_ASSERT( chkVectorS2( a3, 2*11,2*16 ) );


  a2 = a1;

  a3 = (m1 + m2)*(a1 + a2);
  UT_ASSERT( chkVectorS2( a3, 44,64 ) );

  a3 += (m1 + m2)*(a1 + a2);
  UT_ASSERT( chkVectorS2( a3, 2*44,2*64 ) );

  a3 -= (m1 + m2)*(a1 + a2);
  UT_ASSERT( chkVectorS2( a3, 44,64 ) );

  a3 = +((m1 + m2)*(a1 + a2));
  UT_ASSERT( chkVectorS2( a3, 44,64 ) );

//  a3 = (a1 + a2)*(m1 + m2);
//  UT_ASSERT( chkVectorS2( a3, 52,64 ) );

  a3 = (m1 + m2)*(a1 - a2);
  UT_ASSERT( chkVectorS2( a3, 0,0 ) );

//  a3 = (a1 - a2)*(m1 + m2);
//  UT_ASSERT( chkVectorS2( a3, 0,0 ) );

  a3 = (m1 - m2)*(a1 + a2);
  UT_ASSERT( chkVectorS2( a3, 0,0 ) );

//  a3 = (a1 + a2)*(m1 - m2);
//  UT_ASSERT( chkVectorS2( a3, 0,0 ) );

  a3 = (m1 - m2)*(a1 - a2);
  UT_ASSERT( chkVectorS2( a3, 0,0 ) );

//  a3 = (a1 - a2)*(m1 - m2);
//  UT_ASSERT( chkVectorS2( a3, 0,0 ) );
//  UT_ASSERT( chkVectorS2( a1, 1,2 ) );
//  UT_ASSERT( chkVectorS2( a2, 1,2 ) );
//  UT_ASSERT( chkMatrixSymS22( m1, 3,4,5,6 ) );
//  UT_ASSERT( chkMatrixSymS22( m2, 3,4,5,6 ) );

}

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( IO )
{
  //Set the 2nd argument to false to regenerate the pattern file
  output_test_stream output( "IO/dense/MatrixSymS_pattern.txt", true );

  MatrixSymS2 m = Identity();

  output << m << std::endl;
  UT_ASSERT( output.match_pattern() );
  m.dump( 2, output );
  UT_ASSERT( output.match_pattern() );
}

#endif

//############################################################################//
UT_TEST_SUITE_END(MatrixSymD_tests)
