#include "adaptation/metric.h"

#include "geometry/egads/object.h"

#include "library/metric.h"

#include <egads.h>

namespace avro
{

namespace library
{

template<typename type>
MetricField_UniformGeometry<type>::MetricField_UniformGeometry( real_t h , const Field<type,Metric>& field ) :
  MetricField_Analytic( field.topology().number() ),
  FieldInterpolation<type,Metric>(field),
  uniform_( field.topology().number() , h )
{}

template<typename type>
int
MetricField_UniformGeometry<type>::eval( const Points& points , index_t p , const std::vector<index_t>& guesses , Metric& mp )
{
  numerics::SymMatrixD<real_t> m = uniform_(points[p]);
  Entity* entity = points.entity(p);
  avro_assert_msg( entity!=nullptr , "entity is null for vertex %lu?" , p );

  #if 0
  if (entity->number()==2)
  {
    real_t area;
    EGADS::Object* eg = (EGADS::Object*) entity;
    EGADS_ENSURE_SUCCESS( EG_getArea( *eg->object() , NULL , &area ) );

    index_t nelem = area/( std::sqrt(3.)/4.0 );
    real_t h = std::sqrt(area);
    m(0,0) = 1./(h*h);
    m(1,1) = 1./(h*h);
  }
  #endif

  mp.set(m);


  return true;
}

template class MetricField_UniformGeometry<Simplex>;

} // library

} // avro
