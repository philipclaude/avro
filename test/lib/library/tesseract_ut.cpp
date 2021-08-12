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

#include "library/ckf.h"
#include "library/tesseract.h"

#include "graphics/application.h"

using namespace avro;

UT_TEST_SUITE( tesseract_suite )

UT_TEST_CASE( test1 )
{
  std::vector<real_t> x0(4,0);
  std::vector<real_t> length(4,1);
  library::Tesseract tesseract(x0,length);

  tesseract.print();
}
UT_TEST_CASE_END( test1 )

UT_TEST_CASE( test2 )
{
  std::vector<index_t> dims(4,4);
  CKF_Triangulation topology(dims);

  std::vector<real_t> x0(4,0.5);
  std::vector<real_t> length(4,1);
  library::Tesseract tesseract(x0,length);

  topology.points().attach(tesseract);

  topology.points().print(true);
}
UT_TEST_CASE_END( test2 )

void
expand_in_time( const real_t* x , real_t* y ) {

  std::vector<real_t> X(x,x+4);

  // translate to the center of the domain (0.5)^4
  for (coord_t d = 0; d < 4; d++)
    X[d] -= 0.5;

  // scale in time
  real_t t = x[3]; // time between 0 and 1
  real_t s = 1 + t;
  for (coord_t d = 0; d < 3; d++)
    X[d] *= s;
  X[3] *= 5;

  // translate back
  for (coord_t d = 0; d < 4; d++)
    y[d] = X[d] + 0.5;
}

void
squeeze_in_time( const real_t* x , real_t* y ) {

  std::vector<real_t> X(x,x+4);

  // translate to the center of the domain (0.5)^4
  //for (coord_t d = 0; d < 4; d++)
  //  X[d] -= 0.5;

  // scale in time
  real_t t = x[3]; // time between 0 and 1
  real_t s = 1 - 0.5*t;
  X[2] *= s;
  X[3] *= 2;

  // translate back
  for (coord_t d = 0; d < 4; d++)
    y[d] = X[d];// + 0.5;
}


UT_TEST_CASE( test3 )
{
  std::vector<index_t> dims(4,4);
  CKF_Triangulation topology(dims);

  std::vector<real_t> x0(4,0.5);
  std::vector<real_t> length(4,1);
  library::Tesseract tesseract(x0,length);

  tesseract.map_to( &squeeze_in_time , &topology.points() );
  topology.points().attach(tesseract);

  graphics::Viewer app;
  app.add(topology);
  app.run(AVRO_FULL_UNIT_TEST);
}
UT_TEST_CASE_END( test3 )

UT_TEST_SUITE_END( tesseract_suite )
