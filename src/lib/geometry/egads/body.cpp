#include "geometry/egads/body.h"
#include "geometry/egads/object.h"

#include <egads.h>

namespace avro
{

namespace EGADS
{

void
Body::build_hierarchy()
{
  avro_assert( nb_entities()==0 );

  // get the topology of the associated ego with the children
  EGADS_ENSURE_SUCCESS( EG_getTopology( *object_ , &data_.reference , &data_.object_class ,
                        &data_.member_type , data_.data , &data_.nb_children , &data_.children , &data_.senses ) );

  printf("nb_children = %d\n",data_.nb_children);
  // loop through the children obtained from egads
  // and create the topological entities
  for (index_t k=0;k<data_.nb_children;k++)
  {
    // create the new entity
    std::shared_ptr<EGADS::Object> entity;
    entity = std::make_shared<EGADS::Object>(&data_.children[k],this);
    add(entity);

    // build the entity hierarchy
    entity->build_hierarchy();
    entity->set_parent(NULL);
    entity->set_body(this);
  }

  // determine the topological number from the children
  number_ = EGADS::utilities::topological_number(data_.object_class,data_.member_type);

  // assign the parent hierarchy
  build_parents();
}

Body::Entity_ptr
Body::child(index_t k)
{
  avro_assert( k < nb_entities() );
  return entity_[k];
}

void
Body::add_child( ego object , Entity_ptr entity )
{
  children_.insert( {object,entity} );
}

Body::Entity_ptr
Body::lookup( ego object ) const
{
  if (children_.find(object)==children_.end())
    return nullptr;
  return children_.at(object);
}

void
Body::print() const
{
  printf("EGADS: number = %u , class = %s, type = %s at %p\n",number_,
  EGADS::utilities::object_class_name(data_.object_class).c_str(),
  EGADS::utilities::member_type_name(data_.object_class,data_.member_type).c_str(),
  (void*)(this) );

  for (index_t k=0;k<nb_entities();k++)
  {
    printf("\t");
    entity_[k]->print();
  }
}

} // EGADS

} // avro
