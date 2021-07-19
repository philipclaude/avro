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

#include "graphics/application.h"

#include "library/ckf.h"

#include "mesh/field.h"
#include "mesh/field.hpp"
#include "mesh/interpolation.h"

#include "element/quadrature.h"
#include "element/simplex.h"

using namespace avro;

class SomeFunction
{
public:
  SomeFunction()
  {}

  real_t operator() ( const real_t* x ) const
  {
    return 1 + x[0]*x[0] + x[1];
  }
};

UT_TEST_SUITE( interpolation_test_suite )

UT_TEST_CASE( test1 )
{
  CKF_Triangulation topology( {3,3} );
  topology.element().set_basis( BasisFunctionCategory_Lagrange );

  topology.close();
  topology.orient();
  topology.neighbours().compute();
  topology.inverse().build();

  Field<Simplex,real_t> u(topology,2,CONTINUOUS);
  u.build();
  u.element().set_basis( BasisFunctionCategory_Lagrange );

  SomeFunction fcn;
  u.evaluate(fcn);

  //u.print();
  u.dof().print();

  std::shared_ptr< FieldInterpolation<Simplex,real_t> > interpolation;

  interpolation = std::make_shared< FieldInterpolation<Simplex,real_t> >(&u);

  Points points( topology.points().dim() );
  std::vector<real_t> x = {-0.1,0.4};
  points.create(x.data());

  real_t tp;
  interpolation->eval( points , 0 , {0} , tp );

  printf("tp = %g, f = %g \n",tp , fcn(points[0]));

  graphics::Viewer vis;

  //vis.add(topology);
  vis.add(topology);

  vis.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test1 )

UT_TEST_SUITE_END( interpolation_test_suite )
