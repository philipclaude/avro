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
  numerics::SymMatrixD<real_t> m(dim_);// = uniform_(points[p]);
  Entity* entity = points.entity(p);
  avro_assert_msg( entity!=nullptr , "entity is null for vertex %lu?" , p );

  real_t area,h;
  real_t range[4];
  int periodic;
  Entity* face = nullptr;
  if (entity->number()==2)
  {
    face = entity;
  }
  else
  {
    for (index_t k=0;k<entity->nb_parents();k++)
    {
      face = entity->parents(k);
      if (face->number()==2 && face->tessellatable()) break;
    }
    avro_assert(face!=nullptr);
  }

  EGADS::Object* eg = (EGADS::Object*) face;
  EGADS_ENSURE_SUCCESS( EG_getArea( *eg->object() , NULL , &area ) );
  EGADS_ENSURE_SUCCESS( EG_getRange( *eg->object() , range , &periodic ) );

  real_t lu = range[1] - range[0];
  real_t lv = range[3] - range[2];
  real_t area_u = lu*lv;
  h = std::sqrt(area_u);

  h = h*hfactor_;

  m(0,0) = 1./(h*h);
  m(1,1) = 1./(h*h);

  mp.set(m);

  return true;
}

template<typename type>
numerics::SymMatrixD<real_t>
MetricField_UniformGeometry<type>::operator()( const Points& points , index_t p )
{
  Metric m(dim_);
  eval( points , p , {} , m );
  return m;
}

template class MetricField_UniformGeometry<Simplex>;

} // library

} // avro
