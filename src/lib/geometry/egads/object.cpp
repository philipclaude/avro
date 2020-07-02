//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/tools.h"

#include "geometry/egads/object.h"

#include "numerics/coordinate.h"

#include <egads.h>

namespace avro
{

namespace EGADS
{

Object::Object( const Context& context , ego* object ) :
  Entity(0),
  body_(nullptr),
  context_(context),
  object_(object)
{
  std::fill( data_.data , data_.data+4 , 0.0 );
  EGADS_ENSURE_SUCCESS( EG_getInfo( *object_ , &data_.object_class , &data_.member_type ,
                                    &data_.reference , &data_.previous , &data_.next ) );
  number_ = EGADS::utilities::topological_number(data_.object_class,data_.member_type);
  tessellatable_ = EGADS::utilities::object_tessellatable(data_.object_class,data_.member_type);
  interior_ = false;
  sense_required_ = false;
  egads_ = true;
  identifier_ = 0;
  if (data_.member_type == SREVERSE) sign_ = -1;
  name_ = utilities::object_class_name( object_class() ) + "-" + utilities::member_type_name(object_class(),member_type());
}

Object::Object( ego* object , EGADS::Body* body ) :
  Entity(0),
  body_(body),
  context_(body->context()),
  object_(object)
{
  std::fill( data_.data , data_.data+4 , 0.0 );
  EGADS_ENSURE_SUCCESS( EG_getInfo( *object_ , &data_.object_class , &data_.member_type ,
                                    &data_.reference , &data_.previous , &data_.next ) );
  number_ = EGADS::utilities::topological_number(data_.object_class,data_.member_type);
  tessellatable_ = EGADS::utilities::object_tessellatable(data_.object_class,data_.member_type);
  interior_ = false;
  sense_required_ = false;
  egads_ = true;

  identifier_ = EG_indexBodyTopo( *body_->object() , *object_ );

  if (data_.member_type == SREVERSE) sign_ = -1;

  name_ = utilities::object_class_name( object_class() ) + "-" + utilities::member_type_name(object_class(),member_type()) + "-" + std::to_string(identifier_);
}

Object::Object( const Context& context ) :
  Entity(0),
  body_(nullptr),
  context_(context),
  object_(nullptr)
{}

void
Object::construct( ego* object )
{
  set_object(object);
  std::fill( data_.data , data_.data+4 , 0.0 );
  EGADS_ENSURE_SUCCESS( EG_getInfo( *object_ , &data_.object_class , &data_.member_type ,
                                    &data_.reference , &data_.previous , &data_.next ) );
  number_ = EGADS::utilities::topological_number(data_.object_class,data_.member_type);
  tessellatable_ = EGADS::utilities::object_tessellatable(data_.object_class,data_.member_type);
  interior_ = false;
  sense_required_ = false;
  egads_ = true;
  identifier_ = 0;
  if (data_.member_type == SREVERSE) sign_ = -1;
  name_ = utilities::object_class_name( object_class() ) + "-" + utilities::member_type_name(object_class(),member_type());
}

void
Object::set_object( ego* object )
{
  object_ = object;
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

ego
Object::egchild( index_t k ) const
{
  return data_.children[k];
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
Object::inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const
{
  if (number_==0)
  {
    std::fill( u.begin() , u.end() , 0.0 );
    u[0] = 1;
    x[0] = data_.data[0];
    x[1] = data_.data[1];
    if (x.size()==3)
      x[2] = data_.data[2];
    return;
  }

  if (x.size()==2) x.push_back(0);
  avro_assert( x.size()==3 );

  real_t result[3];
  EGADS_ENSURE_SUCCESS( EG_invEvaluate(*object_,x.data(),u.data(),result) );
  for (coord_t d=0;d<3;d++)
    x[d] = result[d];
}

void
Object::inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const
{
  if (number_==0)
  {
    std::fill( u.begin() , u.end() , 0.0 );
    u[0] = 1;
    x[0] = data_.data[0];
    x[1] = data_.data[1];
    x[2] = data_.data[2];
    return;
  }

  if (x.size()==2) x.push_back(0);
  avro_assert( x.size()==3 );

  real_t result[3];
  EGADS_ENSURE_SUCCESS( EG_invEvaluateGuess(*object_,x.data(),u.data(),result) );
  for (coord_t d=0;d<3;d++)
    x[d] = result[d];
}

void
Object::evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const
{
  if (number_==0)
  {
    x[0] = data_.data[0];
    x[1] = data_.data[1];
    x[2] = data_.data[2];
    return;
  }

  real_t data[18];
  EGADS_ENSURE_SUCCESS( EG_evaluate( *object_ , u.data() , data ) );

  for (coord_t d=0;d<3;d++)
    x[d] = data[d];
}

void
Object::print(bool with_children) const
{
  if (with_children)
    for (coord_t i=0;i<body_->number()-number_;i++)
      printf("\t");
  printf("EGADS: number = %u , class = %s, type = %s at %p\n",number_,
  EGADS::utilities::object_class_name(data_.object_class).c_str(),
  EGADS::utilities::member_type_name(data_.object_class,data_.member_type).c_str(),
  (void*)(this) );
  printf("--> data = (%g,%g,%g,%g)\n",data_.data[0],data_.data[1],data_.data[2],data_.data[3]);
  if (!with_children) return;
  for (index_t k=0;k<nb_children();k++)
    child(k).print(with_children);
}

} // EGADS

} // avro
