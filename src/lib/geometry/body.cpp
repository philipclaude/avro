#include "geometry/body.h"
#include "geometry/entity.h"

namespace luna
{

Body::Body( coord_t number ) :
  number_(number)
{}

void
Body::add( Entity_ptr entity )
{
  entity_.push_back(entity);
}

void
Body::build_parents()
{
  for (index_t k=0;k<nb_entities();k++)
    entity_[k]->build_parents();
}

} // luna
