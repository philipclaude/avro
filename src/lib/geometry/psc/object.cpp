#include "geometry/psc/object.h"

namespace avro
{

namespace PSC
{

Object::Object( coord_t number , coord_t dim ) :
  Entity(number),
  dim_(dim)
{}

void
Object::inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const
{
  x[0] = 0;
}

void
Object::inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const
{
  x[0] = 0;
}

void
Object::evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const
{
  x[0] = 0;
}

void
Object::build_hierarchy()
{
  avro_implement;
}

void
Object::print( bool with_children ) const
{
  avro_implement;
}

} // PSC

} // avro
