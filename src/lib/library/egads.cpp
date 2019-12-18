#include "geometry/egads/context.h"
#include "library/egads.h"

namespace luna
{

namespace EGADS
{

Cube::Cube( const Context* context , coord_t number ) :
  Body(*context,&object_)
{}

Cube::Cube( const Context* context , const std::vector<real_t>& lens , const real_t* x0 ) :
  Cube(context,lens.size())
{
  if (lens.size()==2)
  {
    real_t data[6] = { x0[0] , x0[1] , x0[2] , lens[0] , lens[1] , lens[2] };
    EGADS_ENSURE_SUCCESS( EG_makeSolidBody( *context_.get() , BOX , data , &object_ ) );
  }
  else if (lens.size()==3)
  {
    real_t data[6] = { x0[0] , x0[1] , x0[2] , lens[0] , lens[1] , lens[2] };
    EGADS_ENSURE_SUCCESS( EG_makeSolidBody( *context_.get() , BOX , data , &object_ ) );
  }
  this->build_hierarchy();
}

} // EGADS

} // luna
