//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "library/ckf.h"

#include "element/quadrature.h"

#include "mesh/field.hpp"

#include "numerics/integration.h"

using namespace avro;

UT_TEST_SUITE( integration_test_suite )

class Integrand_Monomial : public Integrand<Integrand_Monomial>
{
public:
  typedef real_t T;

public:
  Integrand_Monomial()
  {}

  T operator()( index_t , const QuadraturePoint& , const real_t* x ) const
  {
    return x[0]*x[1];
  }
};

template<typename T>
class Functor_Solution
{
public:
  bool needs_solution() const { return true; }
  bool needs_gradient() const { return false; }
  bool needs_hessian() const { return false; }

  T operator()( const real_t* x ,
                const std::vector<T>& u ,
                const std::vector<T>& ux ,
                const std::vector<T>& uxx ) const
  {
    return u[0];
  }
};

class SomeFunction
{
public:
  SomeFunction( coord_t dim ) :
    dim_(dim)
  {}

  real_t operator() ( const real_t* x ) const
  {
    real_t result = 1;
    for (coord_t d = 0; d < dim_; d++)
      result *= std::exp(x[d])*sin(x[d]);
      //result *= x[d]*sin(x[d]);
    return result;
  }

  real_t analytic() const
  {
    real_t result = 1;
    for (coord_t d = 0; d < dim_; d++)
      result *= 0.5*(1 - std::exp(1)*cos(1) + std::exp(1)*sin(1));
      //result *= sin(1)-cos(1);
    return result;
  }

private:
  coord_t dim_;
};

UT_TEST_CASE( test1 )
{
  CKF_Triangulation topology( {3,3,3,3} );
  topology.element().set_basis( BasisFunctionCategory_Lagrange );

  Integrand_Monomial integrand;

  Functional<Integrand_Monomial> functional(integrand);

  functional.integrate( topology );

  real_t value = functional.value();
  UT_ASSERT_NEAR( value , 0.25 , 1e-12 );

  Field<Simplex,real_t> u(topology,2,CONTINUOUS);
  u.build();
  u.element().set_basis( BasisFunctionCategory_Lagrange );

  SomeFunction fcn(2);
  u.evaluate(fcn);

  typedef Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> Integrand_t;
  Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> integrand2(u);

  Functional<Integrand_t> f2(integrand2);
  f2.integrate( u );
  printf("value = %g\n",f2.value());

}
UT_TEST_CASE_END( test1 )

UT_TEST_CASE( test2 )
{
  std::vector<std::vector<real_t>> error(4,std::vector<real_t>());
  std::vector<std::vector<real_t>> hsize(4,std::vector<real_t>());
  for (coord_t p = 1; p <= 4; p++)
  for (index_t n = 5; n <= 30; n += 5)
  {
    CKF_Triangulation topology( {n,n} );
    topology.element().set_basis( BasisFunctionCategory_Lagrange );

    Field<Simplex,real_t> u(topology,p,CONTINUOUS);
    u.build();
    u.element().set_basis( BasisFunctionCategory_Legendre );

    SomeFunction fcn(2);
    u.evaluate(fcn);

    typedef Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> Integrand_t;
    Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> integrand(u);

    Functional<Integrand_t> f(integrand);
    f.integrate( u );

    real_t analytic = fcn.analytic();
    real_t error0 = analytic - f.value();
    real_t h0 = std::sqrt(1./u.nb_data());
    error[p-1].push_back(error0);
    hsize[p-1].push_back(h0);
    printf("order = %3u, dof = %6lu, h ~ %1.6e, error = %3.6e\n",p, u.nb_data(),h0,error0);
  }

  for (index_t k = 0; k < error.size(); k++)
  {
    index_t order = k +1;
    index_t n0 = error[k].size() -1;
    index_t n1 = error[k].size() -2;

    real_t slope = std::log(error[k][n0]/error[k][n1])/std::log(hsize[k][n0]/hsize[k][n1]);
    printf("slope (p = %3lu) = %1.2f\n",order,slope);

    UT_ASSERT( slope > real_t(order+1) );
  }
}
UT_TEST_CASE_END( test2 )

UT_TEST_CASE( test2_3d )
{
  std::vector<std::vector<real_t>> error(4,std::vector<real_t>());
  std::vector<std::vector<real_t>> hsize(4,std::vector<real_t>());
  for (coord_t p = 1; p <= 4; p++)
  for (index_t n = 2; n <= 8; n+=2)
  {
    CKF_Triangulation topology( {n,n,n} );

    topology.element().set_basis( BasisFunctionCategory_Lagrange );

    Field<Simplex,real_t> u(topology,p,CONTINUOUS);
    u.build();
    u.element().set_basis( BasisFunctionCategory_Bernstein );

    SomeFunction fcn(3);
    u.evaluate(fcn);

    typedef Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> Integrand_t;
    Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> integrand(u);

    Functional<Integrand_t> f(integrand);
    f.integrate( u );

    real_t analytic = fcn.analytic();
    real_t error0 = analytic - f.value();
    real_t h0 = std::pow(u.nb_data(),-1./topology.number());
    error[p-1].push_back(error0);
    hsize[p-1].push_back(h0);
    printf("order = %3u, dof = %6lu, h ~ %1.6e, error = %3.6e\n",p, u.nb_data(),h0,error0);
  }

  for (index_t k=0;k<error.size();k++)
  {
    index_t order = k +1;
    index_t n0 = error[k].size()-1;
    index_t n1 = error[k].size()-2;

    real_t slope = std::log(error[k][n0]/error[k][n1])/std::log(hsize[k][n0]/hsize[k][n1]);
    printf("slope (p = %3lu) = %1.2f\n",order,slope);

    UT_ASSERT( slope > real_t(order+1) );
  }
}
UT_TEST_CASE_END( test2_3d)


UT_TEST_CASE( test2_4d )
{
  std::vector<std::vector<real_t>> error(3,std::vector<real_t>());
  std::vector<std::vector<real_t>> hsize(3,std::vector<real_t>());
  for (coord_t p = 1; p <= 3; p++)
  for (index_t n = 2; n <= 5; n += 1 )
  {
    CKF_Triangulation topology( {n,n,n,n} );
    topology.element().set_basis( BasisFunctionCategory_Lagrange );

    Field<Simplex,real_t> u(topology,p,CONTINUOUS);
    u.build();
    u.element().set_basis( BasisFunctionCategory_Legendre );

    SomeFunction fcn(3);
    u.evaluate(fcn);

    typedef Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> Integrand_t;
    Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> integrand(u);

    Functional<Integrand_t> f(integrand);
    f.integrate( u );

    real_t analytic = fcn.analytic();
    real_t error0 = analytic - f.value();
    real_t h0 = std::pow(u.nb_data(),-1./topology.number());
    error[p-1].push_back(error0);
    hsize[p-1].push_back(h0);
    printf("order = %3u, dof = %6lu, h ~ %1.6e, error = %3.6e\n",p, u.nb_data(),h0,error0);
  }

  for (index_t k=0;k<error.size();k++)
  {
    index_t order = k +1;
    index_t n0 = error[k].size()-1;
    index_t n1 = error[k].size()-2;

    real_t slope = std::log(error[k][n0]/error[k][n1])/std::log(hsize[k][n0]/hsize[k][n1]);
    printf("slope (p = %3lu) = %1.2f\n",order,slope);

    UT_ASSERT( slope > real_t(order+1) );
  }
}
UT_TEST_CASE_END( test2_4d)

UT_TEST_SUITE_END( integration_test_suite )
