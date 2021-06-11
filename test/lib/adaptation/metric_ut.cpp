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

#include "adaptation/metric.h"

#include "common/tools.h"

#include "library/ckf.h"
#include "library/metric.h"

#include "numerics/linear_algebra.h"

using namespace avro;

UT_TEST_SUITE( metric_suite )

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

  numerics::SymMatrixD<real_t> expm = numerics::expm(m);
  UT_ASSERT_NEAR( expm(0,0) , std::exp(1.) , tol );
  UT_ASSERT_NEAR( expm(0,1) , 0. , tol );
  UT_ASSERT_NEAR( expm(1,1) , std::exp(1.) , tol );

  numerics::SymMatrixD<real_t> logm = numerics::logm(m);
  UT_ASSERT_EQUALS( logm(0,0) , 0. );
  UT_ASSERT_EQUALS( logm(0,1) , 0. );
  UT_ASSERT_EQUALS( logm(1,1) , 0. );

  numerics::MatrixD<real_t> q(n,n);
  numerics::VectorD<real_t> lambda(n);
  numerics::eig(m,lambda,q);

  numerics::MatrixD<real_t> m0 = q*diag(lambda)*transpose(q);
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

  numerics::MatrixD<real_t> q(n,n);
  numerics::VectorD<real_t> lambda(n);
  numerics::eig(m,lambda,q);

  numerics::SymMatrixD<real_t> m0 = q*diag(lambda)*transpose(q);
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

  numerics::eig(m,lambda,q);
  m0 = q*diag(lambda)*transpose(q);

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

UT_TEST_CASE(tests_4d)
{
  real_t tol = 1e-12;
  coord_t n = 4;

  Metric m(n);

  // test1: this has a pretty small eigenvalue
  m(0,0) = 2.27;
  m(0,1) = 1.05;
  m(0,2) = 2.45;
  m(0,3) = 1.34;
  m(1,1) = 0.71;
  m(1,2) = 1.15;
  m(1,3) = 0.88;
  m(2,2) = 2.66;
  m(2,3) = 1.56;
  m(3,3) = 1.71;

  numerics::MatrixD<real_t> q(n,n);
  numerics::VectorD<real_t> lambda(n);
  numerics::eig(m,lambda,q);

  numerics::SymMatrixD<real_t> m0 = q*diag(lambda)*transpose(q);
  m0.dump();

  for (index_t i=0;i<n;i++)
  for (index_t j=0;j<n;j++)
    UT_ASSERT_NEAR( m(i,j) , m0(i,j) , tol );

}
UT_TEST_CASE_END(tests_4d)

UT_TEST_CASE(intersect)
{
  for (coord_t d=2;d<=4;d++)
  {
    Metric T(d);

    // test the smr technique for intersection
    Metric T1(d),T2(d);

    T1.set( random_tensor(d) );
    T2.set( random_tensor(d) );

    // intersect with itself
    //T1.intersect(T1,T);
  }
}
UT_TEST_CASE_END(intersect)

UT_TEST_CASE( metric_field_tests )
{
  for (coord_t dim=2;dim<=2;dim++)
  {

    library::MetricField_Uniform function(dim,0.1);

    for (index_t n=4;n<=5;n++)
    {
      std::vector<index_t> sizes(dim,n);
      CKF_Triangulation topology(sizes);
      topology.orient();

      MetricAttachment field(function,topology.points());
      MetricField<Simplex> metric_field(topology,field);

      UT_ASSERT_EQUALS( metric_field.nb_data() , topology.points().nb() );

      metric_field.attachment().set_cells( topology );
      UT_ASSERT( metric_field.check_cells() );

      topology.points().remove( 0 );
      metric_field.remove( 0 );

      UT_ASSERT_EQUALS( metric_field.attachment().nb() , topology.points().nb() );

    }
  }
}
UT_TEST_CASE_END( metric_field_tests )

UT_TEST_SUITE_END( metric_suite )
