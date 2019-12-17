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
    return x[0]*x[1]*sin(x[0])*sin(x[1]);
  }
};

UT_TEST_CASE( test1 )
{
  CKF_Triangulation topology( {3,3} );
  ConicalProductQuadrature quadrature(topology.points().dim());
  quadrature.define();
  topology.master().load_quadrature(quadrature);

  topology.master().set_basis( BasisFunctionCategory_Lagrange );

  Integrand_Monomial integrand;

  Functional<Integrand_Monomial> functional(integrand);

  functional.integrate( topology );

  real_t value = functional.value();
  printf("value = %g\n",value);

  Field<Simplex,real_t> u(topology,2,CONTINUOUS);
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

UT_TEST_CASE( test2 )
{
  //return;
  std::vector<std::vector<real_t>> error(5,std::vector<real_t>());
  std::vector<std::vector<real_t>> hsize(5,std::vector<real_t>());
  for (coord_t p=1;p<=5;p++)
  for (index_t n=5;n<=40;n+=5)
  {
    CKF_Triangulation topology( {n,n} );
    ConicalProductQuadrature quadrature(topology.points().dim());
    quadrature.define();
    topology.master().load_quadrature(quadrature);

    topology.master().set_basis( BasisFunctionCategory_Lagrange );

    Field<Simplex,real_t> u(topology,p,CONTINUOUS);
    u.master().set_basis( BasisFunctionCategory_Lagrange );
    u.master().load_quadrature(quadrature);

    SomeFunction fcn;
    u.evaluate(fcn);

    typedef Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> Integrand_t;
    Integrand_Field<Simplex,real_t,Functor_Solution<real_t>> integrand(u);

    Functional<Integrand_t> f(integrand);
    f.integrate( u );

    real_t analytic = (sin(1)-cos(1))*(sin(1)-cos(1));//1./6;//0.25;
    real_t error0 = analytic - f.value();
    real_t h0 = std::sqrt(1./u.nb_data());
    error[p-1].push_back(error0);
    hsize[p-1].push_back(h0);
    printf("order = %u, dof = %lu, h ~ %g, error = %e\n",p, u.nb(),h0,error0);
  }

  for (index_t k=0;k<error.size();k++)
  {
    index_t order = k +1;
    index_t n0 = error[k].size()-1;
    index_t n1 = error[k].size()-2;

    real_t slope = std::log(error[k][n0]/error[k][n1])/std::log(hsize[k][n0]/hsize[k][n1]);
    printf("slope (p = %lu) = %g\n",order,slope);
    UT_ASSERT( slope >= order+1 );
  }
}
UT_TEST_CASE_END( test2 )

UT_TEST_SUITE_END( integration_test_suite )
