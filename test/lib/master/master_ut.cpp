#include "unit_tester.hpp"

#include "master/simplex.h"
#include "master/quadrature.h"

using namespace luna;

UT_TEST_SUITE( MasterSuite )

UT_TEST_CASE( simplex_tests )
{
  Simplex<Lagrange> master( 2 , 1 );

  ConicalProductQuadrature quadrature( -1 );
  master.loadQuadrature( quadrature );

  master.eval();

}
UT_TEST_CASE_END( simplex_tests )

UT_TEST_SUITE_END( MasterSuite )
