//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "graphics/application.h"

#include <wsss.h>
#include <wsserver.h>

#include <sstream>
#include <unistd.h>

namespace avro
{

namespace graphics
{


WebVisualizer::WebVisualizer()
{
  __web_server__ = this;
}

void
Application<Web_Interface>::connect_client()
{
  json msg;
  msg["commands"] = {"make_tree"};
  msg["make_tree-data"] = scene_->menu();

  wv_broadcastText( const_cast<char*>(msg.dump().c_str()) );
}

void
Application<Web_Interface>::run()
{
  focus_scenes();
  write();

  if (wv_startServer( 7681 , NULL , NULL , NULL , 0 , manager_.context() ) == 0)
  {
    printf("waiting for client...\n");

    #if AVRO_HEADLESS_GRAPHICS
    return;
    #endif
    while (wv_nClientServer(0)==0) {}
    usleep(500000);

    connect_client();

    //receive("colorbar|10");

    while ( wv_statusServer(0) ) usleep(500000);
  }
}

void
Application<Web_Interface>::receive( const std::string& text ) const
{
  printf("received text %s\n",text.c_str());

  std::vector<std::string> message = split(text,"|");
  json response;
  response["command"] = message[0];
  printf("message length = %lu\n",message.size());
  if (message[0]=="colorbar")
  {
    printf("getting colors!\n");
    // loop through the scence graph and get the color limits
    real_t clim[2];
    clim[0] = std::numeric_limits<real_t>::max();
    clim[1] = std::numeric_limits<real_t>::min();
    scene_->get_color_limits( clim );

    printf("color lims = (%g,%g)\n",clim[0],clim[1]);

    response["colorbar-data"] = clim;
    //send(response.dump());
  }
  if (message[0]=="modgprim")
  {

  }
  if (message[0]=="ls")
  {

  }
  if (message[0]=="cd")
  {

  }
}


} // graphics

} // avro
