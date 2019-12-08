#include "common/error.h"

#include "geometry/egads/body.h"
#include "geometry/egads/context.h"
#include "geometry/egads/data.h"
#include "geometry/egads/model.h"

#include <egads.h>

namespace luna
{

namespace EGADS
{

Model::Model( const Context& context ) :
  luna::Model(0),
  context_(context)
{}

Model::Model( const Context& context , const std::string& filename , bool split  ) :
  luna::Model(0),
  context_(context)
{
  int flag = (split) ? 1 : 2;
  int status = EG_loadModel( *context_.get() , flag , filename.c_str() , &object_ );
  luna_assert( status==EGADS_SUCCESS );

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

} // EGADS

} // luna
