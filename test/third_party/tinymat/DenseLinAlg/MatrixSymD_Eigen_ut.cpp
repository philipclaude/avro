// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// MatrixS_btest
// testing of MatrixS<M,N,T> class

#include "unit_tester.hpp"

#include <tools/SANSException.h>

#include "tinymat/dense/dynamic/MatrixSymD.h"
#include "tinymat/dense/dynamic/VectorD.h"
#include "tinymat/dense/dynamic/Eigen.h"
#include "tinymat/dense/dynamic/MatrixD_Diag.h"

#include "chkMatrixD_btest.h"

#include <iostream>
using namespace std;


//Explicitly instantiate the classes so that coverage information is correct
namespace tinymat
{
namespace DLA
{

}
}

using namespace tinymat::DLA;
using namespace tinymat;

typedef double Real;

//############################################################################//
UT_TEST_SUITE( MatrixSymD_test_suite )


//----------------------------------------------------------------------------//
UT_TEST_CASE( Eigen_2x2 )
{
  Real tol = 1e-12;

  MatrixSymD<Real> A1 = { {1.2},
                            {-0.3, 1.5} };
  MatrixSymD<Real> A2 = { {4.9},
                            {0.7, 2.3} };
  MatrixSymD<Real> A3(2);
  A3 = 0;
  MatrixD<Real> E(2,2);
  MatrixSymD<Real> I(2);
  I = Identity();
  VectorD<Real> L(2);
  VectorD<Real> Lnew(2);

  EigenSystem( A1, L, E );
  EigenValues( A1, Lnew );

  UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
  UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );

  VectorD<Real> Zero(2);

  Zero = (A1 - I*L[0])*E.col(0);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );

  Zero = (A1 - I*L[1])*E.col(1);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );

  E = 0;
  EigenVectors( A1, E );

  Zero = (A1 - I*L[0])*E.col(0);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );

  Zero = (A1 - I*L[1])*E.col(1);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );


  EigenSystem( A2, L, E );
  EigenValues( A2, Lnew );

  UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
  UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );

  Zero = (A2 - I*L[0])*E.col(0);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );

  Zero = (A2 - I*L[1])*E.col(1);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );


  // check for divide by zero
  EigenSystem( A3, L, E );
  EigenValues( A3, Lnew );

  UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
  UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );

  UT_ASSERT_EQUALS( 0, L[0] );
  UT_ASSERT_EQUALS( 0, L[1] );

  Zero = (A3 - I*L[0])*E.col(0);

  UT_ASSERT_EQUALS( 0, Zero[0] );
  UT_ASSERT_EQUALS( 0, Zero[1] );

  Zero = (A3 - I*L[1])*E.col(1);

  UT_ASSERT_EQUALS( 0, Zero[0] );
  UT_ASSERT_EQUALS( 0, Zero[1] );
}
UT_TEST_CASE_END( Eigen_2x2 )

//----------------------------------------------------------------------------//
UT_TEST_CASE( Eigen_2x2_ZeroOffDiag )
{
  Real tol = 1e-12;

  MatrixSymD<Real> A = {{2},
                        {0, 4}};
  MatrixD<Real> E(2,2);
  MatrixSymD<Real> I(2);
  I = Identity();
  VectorD<Real> L(2);
  VectorD<Real> Lnew(2);

  EigenSystem( A, L, E );
  EigenValues( A, Lnew );

  UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
  UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );

  VectorD<Real> Zero(2);

  Zero = (A - I*L[0])*E.col(0);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );

  Zero = (A - I*L[1])*E.col(1);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );
}
UT_TEST_CASE_END( Eigen_2x2_ZeroOffDiag )

//----------------------------------------------------------------------------//
UT_TEST_CASE( Eigen_2x2_NearZeroOffDiag )
{
  Real tol = 1e-12;

  MatrixSymD<Real> A = {{2},
                        {1e-44, 2}};
  MatrixD<Real> E(2,2);
  MatrixSymD<Real> I(2);
  I = Identity();
  VectorD<Real> L(2);
  VectorD<Real> Lnew(2);

  EigenSystem( A, L, E );
  EigenValues( A, Lnew );

  UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
  UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );

  VectorD<Real> Zero(2);

  Zero = (A - I*L[0])*E.col(0);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );

  Zero = (A - I*L[1])*E.col(1);

  UT_ASSERT_SMALL( Zero[0], tol );
  UT_ASSERT_SMALL( Zero[1], tol );

  MatrixD<Real> B = E*diag(L)*Transpose(E);

  UT_ASSERT_CLOSE( A(0,0), B(0,0), tol, tol );
  UT_ASSERT_CLOSE( A(1,0), B(1,0), tol, tol );
  UT_ASSERT_CLOSE( A(1,1), B(1,1), tol, tol );
}
UT_TEST_CASE_END( Eigen_2x2_NearZeroOffDiag )

//----------------------------------------------------------------------------//
UT_TEST_CASE( Eigen_Decomposition_2x2 )
{
  Real tol = 1e-12;

  MatrixSymD<Real> A1 = { {1.2},
                            {-0.3, 1.5} };
  MatrixSymD<Real> A2 = { {4.9},
                            {0.7, 2.3} };

  MatrixD<Real> E(2,2);
  MatrixSymD<Real> I(2);
  I = Identity();
  VectorD<Real> L(2);

  EigenSystem( A1, L, E );

  MatrixSymD<Real> A1_decomp = E*tinymat::DLA::diag(L)*tinymat::Transpose(E);

  UT_ASSERT_CLOSE( A1(0,0), A1_decomp(0,0), tol, tol );
  UT_ASSERT_CLOSE( A1(1,0), A1_decomp(1,0), tol, tol );
  UT_ASSERT_CLOSE( A1(1,1), A1_decomp(1,1), tol, tol );

  EigenSystem( A2, L, E );

  MatrixSymD<Real> A2_decomp = E*tinymat::DLA::diag(L)*tinymat::Transpose(E);

  UT_ASSERT_CLOSE( A2(0,0), A2_decomp(0,0), tol, tol );
  UT_ASSERT_CLOSE( A2(1,0), A2_decomp(1,0), tol, tol );
  UT_ASSERT_CLOSE( A2(1,1), A2_decomp(1,1), tol, tol );

}
UT_TEST_CASE_END( Eigen_Decomposition_2x2 )


//----------------------------------------------------------------------------//
UT_TEST_CASE( Eigen_3x3 )
{
  Real tol = 1e-12;

  {
    MatrixSymD<Real> A = { {1.2},
                             {0.5,  1.7},
                             {-0.3, -0.6, 2.3}};

    MatrixD<Real> E(3,3);
    MatrixSymD<Real> I(3);
    I = Identity();
    VectorD<Real> L(3);
    VectorD<Real> Lnew(3);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );

    VectorD<Real> Zero(3);

    Zero = (A - I*L[0])*E.col(0);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    Zero = (A - I*L[1])*E.col(1);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    Zero = (A - I*L[2])*E.col(2);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    MatrixD<Real> B = E*diag(L)*Transpose(E);
    UT_ASSERT_NEAR( B(0,0) , A(0,0) , tol );
    UT_ASSERT_NEAR( B(0,1) , A(0,1) , tol );
    UT_ASSERT_NEAR( B(0,2) , A(0,2) , tol );
    UT_ASSERT_NEAR( B(1,1) , A(1,1) , tol );
    UT_ASSERT_NEAR( B(1,2) , A(1,2) , tol );
    UT_ASSERT_NEAR( B(2,2) , A(2,2) , tol );

  }
  {
    MatrixSymD<Real> A = { {-0.2},
                             {0.8, 0.5},
                             {1.2,-0.4,-1.3}};

    MatrixD<Real> E(3,3);
    MatrixSymD<Real> I(3);
    I = Identity();
    VectorD<Real> L(3);
    VectorD<Real> Lnew(3);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );

    VectorD<Real> Zero(3);

    Zero = (A - I*L[0])*E.col(0);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    Zero = (A - I*L[1])*E.col(1);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    Zero = (A - I*L[2])*E.col(2);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    MatrixD<Real> B = E*diag(L)*Transpose(E);
    UT_ASSERT_NEAR( B(0,0) , A(0,0) , tol );
    UT_ASSERT_NEAR( B(0,1) , A(0,1) , tol );
    UT_ASSERT_NEAR( B(0,2) , A(0,2) , tol );
    UT_ASSERT_NEAR( B(1,1) , A(1,1) , tol );
    UT_ASSERT_NEAR( B(1,2) , A(1,2) , tol );
    UT_ASSERT_NEAR( B(2,2) , A(2,2) , tol );
  }
  {
    MatrixSymD<Real> A = { {2.5091875248975914},
                             {0.52336724167424731, 1.1814971732398194},
                             {-1.2077798854962147, -0.41884286531237125, 1.9665679232991529} };
    MatrixD<Real> E(3,3);
    MatrixSymD<Real> I(3);
    I = Identity();
    VectorD<Real> L(3);
    VectorD<Real> Lnew(3);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );

    VectorD<Real> Zero(3);

    Zero = (A - I*L[0])*E.col(0);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    Zero = (A - I*L[1])*E.col(1);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    Zero = (A - I*L[2])*E.col(2);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );

    MatrixD<Real> B = E*diag(L)*Transpose(E);
    UT_ASSERT_NEAR( B(0,0) , A(0,0) , tol );
    UT_ASSERT_NEAR( B(0,1) , A(0,1) , tol );
    UT_ASSERT_NEAR( B(0,2) , A(0,2) , tol );
    UT_ASSERT_NEAR( B(1,1) , A(1,1) , tol );
    UT_ASSERT_NEAR( B(1,2) , A(1,2) , tol );
    UT_ASSERT_NEAR( B(2,2) , A(2,2) , tol );
  }

  // check for divide by zeros
  {
    MatrixSymD<Real> A(3);
    A = 0;

    MatrixD<Real> E(3,3);
    MatrixSymD<Real> I(3);
    I = Identity();
    VectorD<Real> L(3);
    VectorD<Real> Lnew(3);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );

    UT_ASSERT_EQUALS( 0, L[0] );
    UT_ASSERT_EQUALS( 0, L[1] );
    UT_ASSERT_EQUALS( 0, L[2] );

    VectorD<Real> Zero(3);

    Zero = (A - I*L[0])*E.col(0);

    UT_ASSERT_EQUALS( 0, Zero[0] );
    UT_ASSERT_EQUALS( 0, Zero[1] );
    UT_ASSERT_EQUALS( 0, Zero[2] );

    Zero = (A - I*L[1])*E.col(1);

    UT_ASSERT_EQUALS( 0, Zero[0] );
    UT_ASSERT_EQUALS( 0, Zero[1] );
    UT_ASSERT_EQUALS( 0, Zero[2] );

    Zero = (A - I*L[2])*E.col(2);

    UT_ASSERT_EQUALS( 0, Zero[0] );
    UT_ASSERT_EQUALS( 0, Zero[1] );
    UT_ASSERT_EQUALS( 0, Zero[2] );

    MatrixD<Real> B = E*diag(L)*Transpose(E);
    UT_ASSERT_NEAR( B(0,0) , A(0,0) , tol );
    UT_ASSERT_NEAR( B(0,1) , A(0,1) , tol );
    UT_ASSERT_NEAR( B(0,2) , A(0,2) , tol );
    UT_ASSERT_NEAR( B(1,1) , A(1,1) , tol );
    UT_ASSERT_NEAR( B(1,2) , A(1,2) , tol );
    UT_ASSERT_NEAR( B(2,2) , A(2,2) , tol );
  }

  {
    // See Kopp_2008_Efficient_numerical_diagonalization_of_hermitian_3x3_matrices.pdf
    MatrixSymD<Real> A = { {1e40},
                             {1e19, 1e20},
                             {1e19,  1e9, 1} };
    MatrixD<Real> E(3,3);
    MatrixSymD<Real> I(3);
    I = Identity();
    VectorD<Real> L(3);
    VectorD<Real> Lnew(3);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );

    // These should all be positive, which only happens with Jacobi iterations
    UT_ASSERT_CLOSE( 1e40, L[0], tol, tol );
    UT_ASSERT_CLOSE( 1e20, L[1], tol, tol );
    UT_ASSERT_CLOSE( 0.98, L[2], tol, 1e-10 );

    UT_ASSERT_NEAR(      1., E(0,0), tol );
    UT_ASSERT_NEAR(  1e-21, E(1,0), tol );
    UT_ASSERT_NEAR(  1e-21, E(2,0), tol );

    UT_ASSERT_NEAR( -1e-21, E(0,1), 2e-9 );
    UT_ASSERT_NEAR(      1., E(1,1), tol );
    UT_ASSERT_NEAR(  1e-11, E(2,1), 2e-9 );

    UT_ASSERT_NEAR( -1e-21, E(0,2), 2e-9 );
    UT_ASSERT_NEAR( -1e-11, E(1,2), 2e-9 );
    UT_ASSERT_NEAR(      1., E(2,2), tol );
  }
}
UT_TEST_CASE_END( Eigen_3x3 )

//----------------------------------------------------------------------------//
UT_TEST_CASE( Eigen_4x4 )
{
  Real tol = 1e-12;

  // test 1
  {
    MatrixSymD<Real> A = { {2.3346},
                           {1.1384, 0.7860},
                           {2.5606, 1.2743, 2.8147},
                           {1.4507, 0.9531, 1.6487, 1.8123} };
    MatrixD<Real> E(4,4);
    MatrixSymD<Real> I(4);
    I = Identity();
    VectorD<Real> L(4);
    VectorD<Real> Lnew(4);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );
    UT_ASSERT_CLOSE( L[3], Lnew[3], tol, tol );

    VectorD<Real> Zero(4);

    Zero = (A - I*L[0])*E.col(0);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[1])*E.col(1);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[2])*E.col(2);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[3])*E.col(3);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    MatrixD<Real> B = E*diag(L)*Transpose(E);
    UT_ASSERT_NEAR( B(0,0) , A(0,0) , tol );
    UT_ASSERT_NEAR( B(0,1) , A(0,1) , tol );
    UT_ASSERT_NEAR( B(0,2) , A(0,2) , tol );
    UT_ASSERT_NEAR( B(0,3) , A(0,3) , tol );
    UT_ASSERT_NEAR( B(1,1) , A(1,1) , tol );
    UT_ASSERT_NEAR( B(1,2) , A(1,2) , tol );
    UT_ASSERT_NEAR( B(1,3) , A(1,3) , tol );
    UT_ASSERT_NEAR( B(2,0) , A(2,0) , tol );
    UT_ASSERT_NEAR( B(2,1) , A(2,1) , tol );
    UT_ASSERT_NEAR( B(2,2) , A(2,2) , tol );
    UT_ASSERT_NEAR( B(2,3) , A(2,3) , tol );
    UT_ASSERT_NEAR( B(3,3) , A(3,3) , tol );


  } // end test 1

  // test 2
  {
    MatrixSymD<Real> A = { {2.5647},
                             {1.8781, 2.0246},
                             {1.9452, 1.4695, 1.7409},
                             {1.0231, 1.0652, 1.1118, 0.9585} };
    MatrixD<Real> E(4,4);
    MatrixSymD<Real> I(4);
    I = Identity();
    VectorD<Real> L(4);
    VectorD<Real> Lnew(4);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );
    UT_ASSERT_CLOSE( L[3], Lnew[3], tol, tol );

    VectorD<Real> Zero(4);

    Zero = (A - I*L[0])*E.col(0);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[1])*E.col(1);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[2])*E.col(2);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[3])*E.col(3);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

  } // end test 2

  // test 3
  {
    MatrixSymD<Real> A = { {49.6576},
                             {-14.3583, 35.8262},
                             {7.9538, 13.5091, 37.1327},
                             {-21.7795,21.2868,2.2306,16.8194} };
    MatrixD<Real> E(4,4);
    MatrixSymD<Real> I(4);
    I = Identity();
    VectorD<Real> L(4);
    VectorD<Real> Lnew(4);

    EigenSystem( A, L, E );
    EigenValues( A, Lnew );

    UT_ASSERT_CLOSE( L[0], Lnew[0], tol, tol );
    UT_ASSERT_CLOSE( L[1], Lnew[1], tol, tol );
    UT_ASSERT_CLOSE( L[2], Lnew[2], tol, tol );
    UT_ASSERT_CLOSE( L[3], Lnew[3], tol, tol );

    VectorD<Real> Zero(4);

    Zero = (A - I*L[0])*E.col(0);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[1])*E.col(1);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[2])*E.col(2);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

    Zero = (A - I*L[3])*E.col(3);

    UT_ASSERT_SMALL( Zero[0], tol );
    UT_ASSERT_SMALL( Zero[1], tol );
    UT_ASSERT_SMALL( Zero[2], tol );
    UT_ASSERT_SMALL( Zero[3], tol );

  } // end test 3
}
UT_TEST_CASE_END( Eigen_4x4 )

//############################################################################//
UT_TEST_SUITE_END(MatrixSymD_test_suite)
