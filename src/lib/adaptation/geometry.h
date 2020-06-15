//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_GEOMETRY_H_
#define avro_LIB_ADAPTATION_GEOMETRY_H_

#include "adaptation/primitive.h"

#include "geometry/egads/object.h"

#include "mesh/points.h"

#include "numerics/geometry.h"

#include <array>
#include <vector>

namespace avro
{

// retrieves the geometric parameter coordinates given a facet on the entity
void geometry_params( Entity* e0 , const Points& points , const index_t* v , const index_t nv , real_t* params );

template<typename type>
inline real_t
get_volume( const Topology<type>& topology , Entity* entity , index_t elem , index_t j , const real_t* p )
{
  index_t nf = topology.number()+1;
  std::vector<const real_t*> xk(nf);
  std::vector<real_t> u(2*topology.nv(elem) );
  real_t sign = 1.0;

  coord_t dim = topology.points().dim();
  if (!topology.shape().parameter())
  {
    // neighbour is not in cavity which means we hit a boundary facet
    // set the coordinates
    for (index_t i=0;i<nf;i++)
      xk[i] = topology.points()[ topology( elem, i ) ];
  }
  else
  {
    avro_assert( entity!=nullptr );

    // get the parameter coordinates along the geometry entity
    geometry_params( entity , topology.points() , topology(elem) , topology.nv(elem) , u.data() );

    // set these coordinates into the actual coordinates
    for (index_t i=0;i<nf;i++)
      xk[i] = &u[2*i];

    dim = 2;
    sign = entity->sign();
  }

  // set the last coordinate to the proposed point
  xk[j] = p;

  // check the orientation
  return sign*numerics::simplex_volume(xk,dim);
}

} // avro

#endif
