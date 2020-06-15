//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"

#include "geometry/egads/body.h"
#include "geometry/egads/context.h"
#include "geometry/egads/data.h"
#include "geometry/egads/model.h"

#include <egads.h>

namespace avro
{

namespace EGADS
{

Model::Model( Context* context ) :
  avro::Model(0),
  context_(context),
  mine_(false)
{}

Model::Model( Context* context , const std::string& filename , bool split  ) :
  avro::Model(0),
  context_(context),
  mine_(false)
{
  int flag = (split) ? 1 : 2;
  int status = EG_loadModel( *context_->get() , flag , filename.c_str() , &object_ );
  avro_assert( status==EGADS_SUCCESS );

  // get all bodies
  egoData data;
  status = EG_getTopology( object_ , &data.reference , &data.object_class ,
                                     &data.member_type , NULL , &data.nb_children ,
                                     &data.children , &data.senses );
  printf("detected %d bodies\n",data.nb_children);

  // define and build the hierarchy for each body
  for (int k=0;k<data.nb_children;k++)
  {
    std::shared_ptr<EGADS::Body> body = std::make_shared<EGADS::Body>( data.children+k , this );

    body->build_hierarchy();
    body_.push_back(body);
  }
  determine_number();
}

Model::Model( const std::string& filename , bool split  ) :
  avro::Model(0)
{
  context_ = new Context;
  mine_    = true;

  int flag = (split) ? 1 : 2;
  int status = EG_loadModel( *context_->get() , flag , filename.c_str() , &object_ );
  avro_assert( status==EGADS_SUCCESS );

  // get all bodies
  egoData data;
  status = EG_getTopology( object_ , &data.reference , &data.object_class ,
                                     &data.member_type , NULL , &data.nb_children ,
                                     &data.children , &data.senses );
  printf("detected %d bodies\n",data.nb_children);

  // define and build the hierarchy for each body
  for (int k=0;k<data.nb_children;k++)
  {
    std::shared_ptr<EGADS::Body> body = std::make_shared<EGADS::Body>( data.children+k , this );

    body->build_hierarchy();
    body_.push_back(body);
  }
  determine_number();
}

Model::Model( coord_t number ) :
  avro::Model(number)
{
  context_ = new Context;
  mine_    = true;
}

Entity*
Model::find_entity( index_t id , int object_class ) const
{
  if (nb_bodies()!=1)
    printf("don't know how to find with %lu bodies\n",nb_bodies());

  const Body* b = dynamic_cast<const EGADS::Body*>(&body(0));
	ego object;
	EGADS_ENSURE_SUCCESS( EG_objectBodyTopo( *b->object() , object_class , id , &object ) );
	return b->lookup(object).get();
}

void
Model::determine_number()
{
  number_ = 0;
  for (index_t k=0;k<nb_bodies();k++)
  {
    if (body_[k]->number()>number_)
      number_ = body_[k]->number();
  }
}

void
Model::print() const
{
  for (index_t k=0;k<nb_bodies();k++)
    body_[k]->print();
}

Model::~Model()
{
  if (mine_) delete context_;
}

} // EGADS

} // avro
