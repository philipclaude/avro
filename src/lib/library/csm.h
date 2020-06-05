//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LIBRARY_CSM_H_
#define avro_LIB_LIBRARY_CSM_H_

#include "geometry/egads/model.h"

extern "C"
{
#include <OpenCSM.h>
}
#include <egads.h>

namespace avro
{

class OpenCSM_Model : public EGADS::Model
{
public:
  OpenCSM_Model( const std::string& filename );

  void import();

private:
  modl_T* modl_;
};

} // avro

#endif
