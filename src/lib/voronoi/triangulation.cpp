// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/array.h"

#include "mesh/delaunay/triangulation.h"
#include "mesh/geometrics.h"

#include "triangle/predicates.h"

#include "graph/neighbours.h"

namespace avro
{

namespace delaunay
{

void
Triangulation::compute()
{

  // get the centroid of the vertices
  std::vector<real> xc( vertices_.dim() , 0. );
  geometrics::centroid( linspace<index_t>(0,vertices_.nb()).data() , vertices_.nb() , vertices_ , xc );

  // find the point which is furthest from the centroid
//  index_t kr = geometrics::furthestPoint( xc , vertices_ );
//  index_t dmax
}

void
Triangulation::initialize()
{
  // we need to find a facet to start off with
  // this can be done by looking for a voronoi vertex equidistant to n
  // vertices



}

bool
Triangulation::insphere( const index_t n , const index_t p ) const
{
  real* xn = vertices_[n];
  real* xp = vertices_[p];

  std::vector<real*> xf( number_ );
  for (coord_t d=0;d<number_;d++)
    xf[d] = vertices_[ facet_[d] ];

  if (number_==2)
  {
    avro_implement;
    return ::orient2d(xf[0],xf[1],xf[2]);
  }
  if (number_==3)
  {
    return ::insphere(xf[0],xf[1],xf[2],xn,xp)>0.;
  }
  avro_implement;
  return false;
}

} // delaunay

} // avro
