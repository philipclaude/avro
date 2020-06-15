//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "library/csm.h"

#include "geometry/egads/body.h"

namespace avro
{

OpenCSM_Model::OpenCSM_Model( const std::string& filename0 ) :
  EGADS::Model(2)
{
  std::vector<char> filename(filename0.begin(),filename0.end());
  filename.push_back('\0');
  int status = ocsmLoad( &*filename.begin() , (void**)&modl_ );
  avro_assert(status==0);
  import();
}

void
OpenCSM_Model::import()
{
  int build_to = 0; // execute all branches
  int built_to;
  int nbody;
  int status;

  ocsmCheck(modl_);
  status = ocsmBuild(modl_,build_to,&built_to,&nbody,NULL);

  printf("model has %d bodies\n",modl_->nbody);
  for (int ibody=1;ibody<=modl_->nbody;ibody++)
  {
    if (modl_->body[ibody].onstack!=1) continue;

    // retrieve the egads body
    std::shared_ptr<EGADS::Body> pbody = std::make_shared<EGADS::Body>(&modl_->body[ibody].ebody,this);

    pbody->build_hierarchy();
    body_.push_back(pbody);
  }
  printf("found %lu bodies on stack\n",body_.size());
  determine_number();

}

} // avro
