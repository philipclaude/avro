#include "geometry/model.h"

namespace avro
{

void
Model::get_entities( std::vector<Entity*>& entities ) const
{
  for (index_t k=0;k<body_.size();k++)
    body_[k]->get_entities(entities);
}

} // avro
