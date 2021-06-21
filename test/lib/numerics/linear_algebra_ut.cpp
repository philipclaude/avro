//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "common/tools.h"

#include "numerics/linear_algebra.h"
#include "numerics/mat.h"
#include "numerics/sym.h"

using namespace avro;

matd<real_t>
random_matrix( index_t n )
{
  matd<real_t> X(n,n);
  for (index_t i=0;i<n;i++)
  for (index_t j=0;j<n;j++)
    X(i,j) = random_within(0.,1.);
  return X;
}

symd<real_t>
random_tensor( index_t n )
{
  matd<real_t> X(n,n);
  matd<real_t> Xt(n,n);
  for (index_t i=0;i<n;i++)
  for (index_t j=0;j<n;j++)
  {
    X(i,j)  = random_within(0.,1.);
    Xt(j,i) = X(i,j);
  }

  symd<real_t> N = Xt*X;

  // add n to diagonal for spd-ness
  for (index_t i=0;i<n;i++)
    N(i,i) += n;
  return N;
}

UT_TEST_SUITE(linear_algebra_tests)

UT_TEST_CASE(inverse_tests)
{
  real_t tol = 1e-10;
  index_t ntests = 10;

  for (index_t n=1;n<=4;n++)
  {
    matd<real_t> I(n,n);
    I.eye();

    for (index_t k=0;k<ntests;k++)
    {
      matd<real_t> A(n,n);
      A = random_matrix(n);

      if (std::fabs(numerics::det(A)<1e-12)) continue;

      matd<real_t> Ainv(n,n);
      Ainv = numerics::inverse(A);

      matd<real_t> B(n,n);
      B = A*Ainv;

      matd<real_t> zero(n,n);
      zero = (A*Ainv - I);

      for (index_t i=0;i<n;i++)
      for (index_t j=0;j<n;j++)
        UT_ASSERT_NEAR( zero(i,j) , 0.0 , tol );

    }
  }

  matd<real_t> A(5,5),Ainv(5,5),B(5,4);
  UT_CATCH_EXCEPTION( Ainv = numerics::inverse(A) );
  UT_CATCH_EXCEPTION( Ainv = numerics::inverse(B) );
}
UT_TEST_CASE_END(inverse_tests)

UT_TEST_CASE( symd_eign_test )
{
  real_t tol = 1e-12;

  for (index_t n = 2; n < 10; n++) {

    symd<real_t> A = random_tensor(n);

    std::pair< vecd<real_t> , matd<real_t> > e = numerics::eign(A);
    const matd<real_t>& Q = e.second;
    const vecd<real_t>& L = e.first;

    matd<real_t> A1 = Q * numerics::diag(L) * numerics::transpose(Q);
    for (index_t i = 0; i < n; i++)
    for (index_t j = 0; j < n; j++)
      UT_ASSERT_NEAR( A(i,j) , A1(i,j) , tol );
  }

}
UT_TEST_CASE_END( symd_eign_test )


UT_TEST_SUITE_END(linear_algebra_tests)
