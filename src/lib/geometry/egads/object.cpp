#include "geometry/egads/object.h"

#include "numerics/coordinate.h"

#include <egads.h>

namespace luna
{

namespace EGADS
{

Object::Object( const Context& context , ego* object ) :
  Entity(0),
  body_(nullptr),
  context_(context),
  object_(object)
{
  EGADS_ENSURE_SUCCESS( EG_getInfo( *object_ , &data_.object_class , &data_.member_type ,
                                    &data_.reference , &data_.previous , &data_.next ) );
  number_ = EGADS::utilities::topological_number(data_.object_class,data_.member_type);
}

Object::Object( ego* object , EGADS::Body* body ) :
  Entity(0),
  body_(body),
  context_(body->context()),
  object_(object)
{
  EGADS_ENSURE_SUCCESS( EG_getInfo( *object_ , &data_.object_class , &data_.member_type ,
                                    &data_.reference , &data_.previous , &data_.next ) );
  number_ = EGADS::utilities::topological_number(data_.object_class,data_.member_type);
}

ego*
Object::object()
{
  return object_;
}

const ego*
Object::object() const
{
  return object_;
}

void
Object::build_hierarchy()
{
  EG_getTopology( *object_ , &data_.reference , &data_.object_class  , &data_.member_type ,
                  data_.data , &data_.nb_children , &data_.children , &data_.senses );

  for (int k=0;k<data_.nb_children;k++)
  {
    // get the index of this child in the body's full list of children
    std::shared_ptr<Entity> entity = body_->lookup( data_.children[k] );

    if (entity == nullptr)
    {
      // the body was not aware of this entity, create a new one and add it to the body
      entity = std::make_shared<EGADS::Object>( &data_.children[k] , body_ );

      // add the new entity to the body's full list
      body_->add_child( data_.children[k] , entity );
    }

    // look to see if the child edge is repeated twice in the loop
    // this means the edge needs a sense for EG_getEdgeUV
    if (data_.object_class==LOOP)
    {
      for (int i=0;i<data_.nb_children;i++)
      {
        if (i==k) continue;
        if (data_.children[i]==data_.children[k])
        {
          entity->set_sense_required(true);
          break;
        }
      }
    }

    // add the new entity as a child of this one
    Tree<Entity>::add_child(entity);

    // build the child
    entity->build_hierarchy();
    entity->set_parent(this);
  }
}

void
Object::inverse( numerics::Coordinate& x , numerics::Coordinate& u ) const
{
  x[0] = 0;
}

void
Object::evaluate( const numerics::Coordinate& u , numerics::Coordinate& x ) const
{
  x[0] = 0;
}

void
Object::print() const
{
  for (index_t i=0;i<body_->number()-number_;i++)
    printf("\t");
  printf("EGADS: number = %u , class = %s, type = %s at %p\n",number_,
  EGADS::utilities::object_class_name(data_.object_class).c_str(),
  EGADS::utilities::member_type_name(data_.object_class,data_.member_type).c_str(),
  (void*)(this) );
  for (index_t k=0;k<nb_children();k++)
    child(k).print();
}

} // EGADS

} // luna
