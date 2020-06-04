#include "unit_tester.hpp"

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/field.h"
#include "mesh/field.hpp"
#include "mesh/interpolation.h"

#include "master/quadrature.h"
#include "master/simplex.h"

using namespace avro;

class SomeFunction
{
public:
  SomeFunction( coord_t dim ) :
    dim_(dim)
  {}

  real_t operator() ( const real_t* x ) const
  {
    return 1 + x[0]*x[0] + x[1];
    //real_t result = 1;
    //for (coord_t d=0;d<dim_;d++)
    //  result *= std::exp(x[d])*sin(x[d]);
    //return result;
  }

private:
  coord_t dim_;
};

UT_TEST_SUITE( interpolation_test_suite )

UT_TEST_CASE( test1 )
{
  CKF_Triangulation topology( {3,3} );
  ConicalProductQuadrature quadrature(topology.points().dim());
  quadrature.define();
  topology.master().load_quadrature(quadrature);
  topology.master().set_basis( BasisFunctionCategory_Lagrange );

  topology.close();
  topology.orient();
  topology.neighbours().compute();
  topology.inverse().build();

  Field<Simplex,real_t> u(topology,2,CONTINUOUS);
  u.build();
  u.master().set_basis( BasisFunctionCategory_Lagrange );
  u.master().load_quadrature(quadrature);

  SomeFunction fcn(2);
  u.evaluate(fcn);

  //u.print();
  u.dof().print();

  std::shared_ptr< FieldInterpolation<Simplex,real_t> > interpolation;

  interpolation = std::make_shared< FieldInterpolation<Simplex,real_t> >(u);

  Points points( topology.points().dim() );
  std::vector<real_t> x = {-0.1,0.4};
  points.create(x.data());

  real_t tp;
  interpolation->eval( points , 0 , {0} , tp );

  printf("tp = %g, f = %g \n",tp , fcn(points[0]));

  graphics::Visualizer vis;

  //vis.add_topology(topology);
  vis.add_topology(topology);

  // test the EPS export
  vis.run();
/*
  interpolation = std::make_shared< GeometryMetric<Simplex> >(u);
  interpolation->eval( topology.points() , 0 , {0} , tp );
*/

}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( interpolation_test_suite )
