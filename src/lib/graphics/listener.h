//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_LISTENER_H_
#define avro_LIB_GRAPHICS_LISTENER_H_

#include "common/directory.h"
#include "common/json.h"

namespace avro
{

namespace graphics
{


class Listener
{
public:
  Listener() {}

  void send( const json& request , json& response );

  std::string pwd() const { return directory_.pwd(); }

private:
  Directory directory_;
};

} // graphics

} // avro

#endif
