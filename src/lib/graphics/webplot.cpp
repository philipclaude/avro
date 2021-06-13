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
Application<Web_Interface>::run( const std::string& view )
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

void
Application<Web_Interface>::receive( const std::string& text )
{
  //printf("received text %s\n",text.c_str());

  std::vector<std::string> message = split(text,"|");

  std::vector<std::string> commands;
  commands.push_back(message[0]);

  json response;
  response["commands"] = commands;

  if (message[0]=="colorbar")
  {
    real_t clim[2];
    clim[0] = std::numeric_limits<real_t>::max();
    clim[1] = std::numeric_limits<real_t>::min();

    // generate color values
    int ncol = atoi(message[1].c_str());
    std::vector<float> values;
    colormap_.generate(ncol,values);

    int root = atoi(message[2].c_str());
    if (root>=0)
      scene_->primitive(root).get_field_limits( clim );
    else
    {
      float lims[2];
      colormap_.get_limits( lims );
      clim[0] = lims[0];
      clim[1] = lims[1];
    }

    printf("color lims = (%g,%g)\n",clim[0],clim[1]);

    json json_data;
    json_data["values"] = values;
    json_data["limits"] = {clim[0],clim[1]};

    response["colorbar-data"] = json_data;
    send(response.dump());
  }
  if (message[0]=="select-field")
  {
    // parse the incoming message
    index_t root = atoi(message[1].c_str());
    std::string id0 = message[3];
    std::vector<std::string> s = split(id0,"-");
    std::string id = s[0];
    index_t rank = atoi(s[1].c_str());

    // assign which field and rank is active
    if (id=="0x0")
      scene_->primitive(root).set_active(id,0);
    else
    {
      std::string name = scene_->primitive(root).topology().fields().id2name(id);
      scene_->primitive(root).set_active(name,rank);
    }

    // re-extract the geometric primitives (mostly just assigning the colors)
    scene_->primitive(root).extract();
    for (index_t k=0;k<scene_->primitive(root).nb_children();k++)
      scene_->primitive(root).child(k).extract();

    // overwrite the primitives to the graphics manager
    write();

    // adjust the colorbar
    receive("colorbar|10|"+std::to_string(root));
  }
  if (message[0]=="plotclip")
  {
    std::string dir0 = message[1];
    std::vector<std::string> data = split(message[2],",");

    int dir = 1;
    if (dir0=="L") dir *= -1;

    if (dir0=="H")
    {
      // request to hide the clip
      clip_plane_.hide(manager_);
      return;
    }

    // update the plane
    real_t d = stod(data[0]);
    real_t angles[2] = {stod(data[1]),stod(data[2])};
    clip_plane_.update(d,angles,dir);
    clip_plane_.plot(manager_,focus_,glm::identity());
  }
  if (message[0]=="clip")
  {
    scene_->write( manager_ , &clip_plane_ );
  }
  if (message[0]=="ls")
  {

  }
  if (message[0]=="cd")
  {

  }
  if (message[0]=="colormap")
  {
    colormap_.change_style(message[1]);
    write();
  }
}


} // graphics

} // avro
