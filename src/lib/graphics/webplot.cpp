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

void
Application<Web_Interface>::connect_client()
{
  json msg;
  msg["commands"] = {"make_tree"};
  msg["make_tree-data"] = scene_.menu();

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

    while ( wv_statusServer(0) ) usleep(500000);
  }
}

} // graphics

} // avro
