#include "unit_tester.hpp"

#include "adaptation/metric.h"

#include "common/tools.h"

#include "numerics/linear_algebra.h"

using namespace luna;

UT_TEST_SUITE( Metric_suite )

numerics::SymMatrixD<real_t>
random_tensor( index_t n )
{
  numerics::MatrixD<real_t> X(n,n);
  numerics::MatrixD<real_t> Xt(n,n);
  for (index_t i=0;i<n;i++)
  for (index_t j=0;j<n;j++)
  {
    X(i,j)  = random_within(0.,1.);
    Xt(j,i) = X(i,j);
  }

  numerics::SymMatrixD<real_t> N = Xt*X;

  // add n to diagonal for spd-ness
  for (index_t i=0;i<n;i++)
    N(i,i) += n;
  return N;
}

numerics::SymMatrixD<real_t>
identity( index_t n )
{
  numerics::SymMatrixD<real_t> m(n);
  for (index_t i=0;i<n;i++)
    m(i,i) = 1;
  return m;
}

UT_TEST_CASE(Metric_tests_2d)
{
  real_t tol = 1e-12;

  coord_t n = 2;

  // 2d
  Metric m(n);
  m.set(identity(n));

  numerics::SymMatrixD<real_t> expm = numerics::exp(m);
  UT_ASSERT_NEAR( expm(0,0) , std::exp(1.) , tol );
  UT_ASSERT_NEAR( expm(0,1) , 0. , tol );
  UT_ASSERT_NEAR( expm(1,1) , std::exp(1.) , tol );

  numerics::SymMatrixD<real_t> logm = numerics::log(m);
  UT_ASSERT_EQUALS( logm(0,0) , 0. );
  UT_ASSERT_EQUALS( logm(0,1) , 0. );
  UT_ASSERT_EQUALS( logm(1,1) , 0. );

  numpack::DLA::MatrixD<real_t> q(n,n);
  numerics::VectorD<real_t> lambda(n);
  numpack::DLA::EigenSystem(m,lambda,q);

  numerics::MatrixD<real_t> m0 = q*numpack::DLA::diag(lambda)*numpack::Transpose(q);
  m0.dump();

  UT_ASSERT_EQUALS( m0.m() , m.m() );
  UT_ASSERT_EQUALS( m0.n() , m.n() );

  for (index_t i=0;i<n;i++)
  for (index_t j=0;j<n;j++)
    UT_ASSERT_NEAR( m(i,j) , m0(i,j) , tol );
}
UT_TEST_CASE_END(Metric_tests_2d)

UT_TEST_CASE(Metric_tests_3d)
{
  real_t tol = 1e-12;
  coord_t n = 3;

  Metric m(n);

  // test1
  m(0,0) = 1.2;
  m(0,1) = 0.5;
  m(0,2) = 1.7;
  m(1,1) = -0.3;
  m(1,2) = -0.6;
  m(2,2) = 2.3;

  numpack::DLA::MatrixD<real_t> q(n,n);
  numerics::VectorD<real_t> lambda(n);
  numpack::DLA::EigenSystem(m,lambda,q);

  numerics::SymMatrixD<real_t> m0 = q*numpack::DLA::diag(lambda)*numpack::Transpose(q);
  m0.dump();

  for (index_t i=0;i<n;i++)
  for (index_t j=0;j<n;j++)
    UT_ASSERT_NEAR( m(i,j) , m0(i,j) , tol );

  // test2
  m(0,0) = -.2;
  m(0,1) = 0.8;
  m(0,2) = 0.5;
  m(1,1) = 1.2;
  m(1,2) = -0.4;
  m(2,2) = -1.3;

  numpack::DLA::EigenSystem(m,lambda,q);

  m0 = q*numpack::DLA::diag(lambda)*numpack::Transpose(q);

  for (index_t i=0;i<n;i++)
  for (index_t j=0;j<n;j++)
    UT_ASSERT_NEAR( m(i,j) , m0(i,j) , tol );

  Metric m1(n),m2(n);
  m1.set( identity(n) );
  m2.set( identity(n) );

  Metric m3(n);
  std::vector<real_t> alpha(2,0.5);
  interp({0.5,0.5},{m1,m2},m3);
  m3.dump();
}
UT_TEST_CASE_END(Metric_tests_3d)

UT_TEST_SUITE_END( Metric_suite )
