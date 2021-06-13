//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/metric.h"

#include "geometry/egads/object.h"

#include "library/metric.h"

#include <egads.h>

namespace avro
{

namespace library
{

template<typename type>
MetricField_UniformGeometry<type>::MetricField_UniformGeometry( coord_t dim , real_t h ) :
  MetricField_Analytic(dim),
  FieldInterpolation<type,Metric>(nullptr),
  hfactor_(h)
{
  analytic_ = true;
}

template<typename type>
int
MetricField_UniformGeometry<type>::eval( const Points& points , index_t p , const std::vector<index_t>& guesses , Metric& mp )
{
  symd<real_t> m(dim_);// = uniform_(points[p]);
  Entity* entity = points.entity(p);
  avro_assert_msg( entity!=nullptr , "entity is null for vertex %lu?" , p );

  if (entity->number()==2)
  {
    // evaluate the metric on the face
    eval_face(points,p,entity,mp);
    return 0;
  }

  // interpolate the metrics evaluated on each face
  std::vector<symd<real_t>> metrics;
  for (index_t k=0;k<entity->nb_parents();k++)
  {
    Entity* face = entity->parents(k);
    if (face->number()==2 && face->tessellatable())
    {
      Metric mk(mp.n());
      eval_face( points , p , face , mk );
      metrics.push_back(mk);
    }
  }
  avro_assert( metrics.size()!=0 );
  std::vector<real_t> alpha(metrics.size(),1./metrics.size());
  mp = numerics::interp( alpha , metrics );
  return 0;
}

template<typename type>
int
MetricField_UniformGeometry<type>::eval_face( const Points& points , index_t p , Entity* entity , Metric& mp )
{
  avro_assert( entity->number()==2 );

  symd<real_t> m(dim_);
  real_t area,h;
  real_t range[4];
  int periodic;

  EGADS::Object* eg = (EGADS::Object*) entity;
  #ifndef AVRO_NO_ESP
  EGADS_ENSURE_SUCCESS( EG_getArea( *eg->object() , NULL , &area ) );
  #else
  area = 1.0;
  #endif
  EGADS_ENSURE_SUCCESS( EG_getRange( *eg->object() , range , &periodic ) );

  real_t lu = range[1] - range[0];
  real_t lv = range[3] - range[2];
  real_t area_u = lu*lv;
  h = std::sqrt(area_u);

  h = h*hfactor_;

  m(0,0) = 1./(h*h);
  m(1,1) = 1./(h*h);
  mp.set(m);

  return 0; // returned element should not be used
}

template<typename type>
symd<real_t>
MetricField_UniformGeometry<type>::operator()( const Points& points , index_t p )
{
  Metric m(dim_);
  if (p < points.nb_ghost())
  {
    m = 0.0;
    return m;
  }
  eval( points , p , {} , m );
  return m;
}

template class MetricField_UniformGeometry<Simplex>;

} // library

} // avro
