//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/listener.h"

#include <json/json.hpp>

namespace avro
{

namespace graphics
{

void
Listener::send( const json& request , json& response )
{
  std::string command = request.at("command");

  if (command=="ls")
  {
    std::vector<json> items;
    directory_.ls(items);
    response["ls-response"] = items;
  }
  if (command=="cd")
  {
    directory_.cd( request["data"] );
  }
}

} // graphics

} // avro
