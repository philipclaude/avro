#include "unit_tester.hpp"

#include "master/quadrature.h"

#include "numerics/functions.h"
#include "numerics/integration.h"

using namespace avro;

UT_TEST_SUITE( quadrature_test_suite )

template<int dim>
struct monomial : public Integrand<monomial<dim>>
{
  void operator() ( const real_t* x , real_t& f ) const
  {
    f = 1.;
    for (coord_t d=0;d<dim;d++)
      f *= pow(x[d],E[d]);
  }

  real_t analytic()
  {
    real_t num = 1.;
    real_t den = 0.;
    for (coord_t d=0;d<dim;d++)
    {
      num *= numerics::factorial(E[d]);
      den += E[d];
    }
    return num/(numerics::factorial(den+dim));
  }

  index_t E[dim];
};

template<int dim>
real_t
evaluate( Quadrature& quadrature , const monomial<dim>& f )
{
  index_t order = 0;
  for (coord_t d=0;d<dim;d++)
    order += f.E[d];
  quadrature.order() = order;
  quadrature.define();

  real_t I = 0.;
  real_t df;
  real_t wsum = 0.;
  for (index_t k=0;k<quadrature.nb();k++)
  {
    f( quadrature.x(k) , df );
    I += quadrature.w(k)*df;
    wsum += quadrature.w(k);

    //if (dim==2)
    //  printf("x[%lu] = (%g,%g)\n",k,quadrature.x(k)[0],quadrature.x(k)[1]);
  }
  if ( fabs(wsum-1.0)>1e-12 ) printf("weights should sum to 1!! wsum = %g\n",wsum);
  I *= quadrature.volume(); // jacobian

  return I;
}

UT_TEST_CASE( line_cp_tests )
{
  const real_t tol = 1e-12;

  ConicalProductQuadrature quadrature( 1 , -1 );

  monomial<1> func;
  real_t I = 0.;

  func.E[0] = 1; // linear
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 2;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 3;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 4;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 5;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 6;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 7;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 8;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 9;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 10;
  I = evaluate<1>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );
}
UT_TEST_CASE_END( line_cp_tests )

UT_TEST_CASE( triangle_cp_tests )
{
  const real_t tol = 1e-12;

  ConicalProductQuadrature quadrature( 2 , -1 );

  monomial<2> func;
  real_t I = 0.;

  func.E[0] = 1;
  func.E[1] = 0;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 0;
  func.E[1] = 1;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 1;
  func.E[1] = 1;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 1;
  func.E[1] = 2;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 2;
  func.E[1] = 1;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 2;
  func.E[1] = 2;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 2;
  func.E[1] = 3;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 3;
  func.E[1] = 2;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 3;
  func.E[1] = 3;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 3;
  func.E[1] = 4;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 4;
  func.E[1] = 3;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 4;
  func.E[1] = 4;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 4;
  func.E[1] = 5;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 5;
  func.E[1] = 4;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 5;
  func.E[1] = 5;
  I = evaluate<2>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );
}
UT_TEST_CASE_END( triangle_cp_tests )

UT_TEST_CASE( tetrahedron_cp_tests )
{
  const real_t tol = 1e-12;

  ConicalProductQuadrature quadrature( 3 , -1 );

  monomial<3> func;
  real_t I = 0.;

  func.E[0] = 1;
  func.E[1] = 0;
  func.E[2] = 0;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 0;
  func.E[1] = 1;
  func.E[2] = 0;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 0;
  func.E[1] = 0;
  func.E[2] = 1;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 1;
  func.E[1] = 1;
  func.E[2] = 0;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 1;
  func.E[1] = 0;
  func.E[2] = 1;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 0;
  func.E[1] = 1;
  func.E[2] = 1;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 1;
  func.E[1] = 1;
  func.E[2] = 1;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 2;
  func.E[1] = 0;
  func.E[2] = 0;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 0;
  func.E[1] = 2;
  func.E[2] = 0;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 0;
  func.E[1] = 0;
  func.E[2] = 2;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 2;
  func.E[1] = 1;
  func.E[2] = 0;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 2;
  func.E[1] = 0;
  func.E[2] = 1;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 1;
  func.E[1] = 2;
  func.E[2] = 0;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  func.E[0] = 0;
  func.E[1] = 2;
  func.E[2] = 1;
  I = evaluate<3>( quadrature , func );
  printf("I = %g, analytic = %g\n",I,func.analytic());
  UT_ASSERT_NEAR( I , func.analytic() , tol );

  for (index_t i=0;i<6;i++)
  for (index_t j=0;j<6;j++)
  for (index_t k=0;k<6;k++)
  {
    if (i==0 && j==0 && k==0) continue;
    func.E[0] = i;
    func.E[1] = j;
    func.E[2] = k;
    I = evaluate<3>( quadrature , func );
    //printf("I = %g, analytic = %g\n",I,func.analytic());
    UT_ASSERT_NEAR( I , func.analytic() , tol );
  }
}
UT_TEST_CASE_END( tetrahedron_cp_tests )

UT_TEST_CASE( pentatope_cp_tests )
{
  const real_t tol = 1e-12;

  ConicalProductQuadrature quadrature( 4 , -1 );

  monomial<4> func;
  real_t I = 0.;

  for (index_t i=0;i<5;i++)
  for (index_t j=0;j<5;j++)
  for (index_t k=0;k<5;k++)
  for (index_t m=0;m<5;m++)
  {
    if (i==0 && j==0 && k==0) continue;
    func.E[0] = i;
    func.E[1] = j;
    func.E[2] = k;
    func.E[3] = m;
    I = evaluate<4>( quadrature , func );
    //printf("I = %g, analytic = %g\n",I,func.analytic());
    UT_ASSERT_NEAR( I , func.analytic() , tol );
  }

}
UT_TEST_CASE_END( pentatope_cp_tests )

UT_TEST_SUITE_END( quadrature_test_suite )
