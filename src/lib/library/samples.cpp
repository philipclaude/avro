//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/tools.h"

#include "library/samples.h"

namespace avro
{

namespace library
{

TwoTriangles::TwoTriangles() :
  Topology<Simplex>(points_,2),
  points_(3),
  edges_(points_,1)
{
  this->set_sorted(false);

  real_t x0[3] = {0,0,0};
  real_t x1[3] = {1,0,0};
  real_t x2[3] = {0,1,0};
  real_t x3[3] = {1,1,0};
  points_.create( x0 );
  points_.create( x1 );
  points_.create( x2 );
  points_.create( x3 );

  index_t t0[3] = {0,1,2};
  index_t t1[3] = {1,3,2};

  add( t0 , 3 );
  add( t1 , 3 );

  index_t e0[2] = {0,1};
  index_t e1[2] = {1,2};
  index_t e2[2] = {2,0};
  index_t e3[2] = {1,3};
  index_t e4[2] = {3,2};

  edges_.add( e0 , 2 );
  edges_.add( e1 , 2 );
  edges_.add( e2 , 2 );
  edges_.add( e3 , 2 );
  edges_.add( e4 , 2 );
}

RegularPolygon::RegularPolygon( index_t nb_side ) :
  Topology<Polytope>(points_,2),
  points_(2)
{
  real_t dtheta = 2.0*M_PI/real_t(nb_side);

  std::vector<real_t> x(2);
  for (index_t k=0;k<nb_side;k++)
  {
    x[0] = cos(k*dtheta);
    x[1] = sin(k*dtheta);

    points_.create(x.data());
  }

  Table<int>& incidence = points_.incidence();
  int hrep[2];
  for (index_t k=0;k<nb_side;k++)
  {
    hrep[0] = k;
    if (k==nb_side-1) hrep[1] = 0;
    else hrep[1] = k+1;
    incidence.add(hrep,2);
  }

  std::vector<index_t> polytope = linspace(nb_side);
  add(polytope.data(),polytope.size());

}

} // library

} // avro
