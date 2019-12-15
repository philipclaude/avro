#include "unit_tester.hpp"

#include "library/ckf.h"

#include "master/quadrature.h"

#include "mesh/field.hpp"

#include "numerics/integration.h"

using namespace luna;

UT_TEST_SUITE( integration_test_suite )

class Integrand_Monomial : public Integrand<Integrand_Monomial>
{
public:
  typedef real_t T;

public:
  Integrand_Monomial()
  {}

  T operator()( index_t , const real_t* , const real_t* x ) const
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
  real_t operator() ( const real_t* x ) const
  {
    return x[0]*x[1];
  }
};

UT_TEST_CASE( test1 )
{
  CKF_Triangulation topology( {100,100} );
  ConicalProductQuadrature quadrature(topology.points().dim(),2);
  quadrature.define();
  topology.master().load_quadrature(quadrature);

  topology.master().set_basis( BasisFunctionCategory_Lagrange );

  Integrand_Monomial integrand;

  Functional<Integrand_Monomial> functional(integrand);

  functional.integrate( topology );

  real_t value = functional.value();
  printf("value = %g\n",value);

  Field<Simplex,real_t> u(topology,3,CONTINUOUS);
  u.master().set_basis( BasisFunctionCategory_Lagrange );
  u.master().load_quadrature(quadrature);

  SomeFunction fcn;
  u.evaluate(fcn);

  typedef Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> Integrand_t;
  Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> integrand2(u);

  Functional<Integrand_t> f2(integrand2);
  f2.integrate( u );

  printf("value = %g\n",f2.value());
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( integration_test_suite )
