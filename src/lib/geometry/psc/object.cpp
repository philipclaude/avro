#include "common/tools.h"

#include "geometry/psc/object.h"

namespace avro
{

namespace PSC
{

Object::Object( coord_t number , coord_t dim ) :
  Entity(number),
  dim_(dim)
{
  tessellatable_ = true;
  interior_ = false;
}

void
Object::build_hierarchy()
{
  // not necesssary but declared virtual...maybe should just
  // use in egads geometry
  avro_assert_not_reached;
}

void
Object::print( bool with_children ) const
{
  printf("PSC entity with number %hu at %p\n",number_,(void*)(this));
  if (!with_children) return;
  for (index_t k=0;k<nb_children();k++)
    child(k).print(with_children);
}

} // PSC

} // avro
