//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GRAPHICS_CLIENT_H_
#define avro_LIB_GRAPHICS_CLIENT_H_

#include "common/error.h"
#include "common/tools.h"
#include "common/types.h"

#include "graphics/clipping.h"
#include "graphics/colormap.h"
#include "graphics/gl.h"
#include "graphics/manager.h"
#include "graphics/scene.h"
#include "graphics/window.h"
#include "graphics/wv.h"

#include <map>
#include <memory>
#include <unordered_set>
#include <string>
#include <vector>

namespace avro
{

class Mesh;
class TopologyBase;

namespace graphics
{

class Widget;

class ApplicationBase
{
public:

  void write();
  Colormap& colormap() { return colormap_; }

protected:
  ApplicationBase( GraphicsManager& manager ) :
    manager_(manager)
  {}

  virtual ~ApplicationBase() {}

  virtual void run() = 0;

  void focus_scenes();

protected:
  std::vector<SceneGraph*> scenes_;
  Colormap colormap_;

private:
  GraphicsManager& manager_;

protected:
  real_t bounding_box_[6]; // (xmin,ymin,zmin,xmax,ymax,zmax)
  real_t focus_[4]; // (xcenter,ycenter,zcenter,scale)
};

template<typename type> class Application;

template<typename T> struct GLFW_Interface;
struct Web_Interface;

#ifdef AVRO_WITH_GL

template<typename API_t>
class Application<GLFW_Interface<API_t>> : public ApplicationBase
{
public:
  Application();

  void initialize();
  void run();

  bool& restart() { return restart_; }

protected:
  void add_window( GLFW_Window* window )
    { window_.push_back(window); }

protected:
  API_t manager_;
  std::vector<GLFW_Window*> window_;
  bool restart_;
};

#endif

template<>
class Application<Web_Interface> : public ApplicationBase
{
public:
  Application() :
    ApplicationBase(manager_),
    clip_plane_(3)
  {
    scene_ = std::make_shared<SceneGraph>();
    scenes_.push_back(scene_.get());
    scene_->set_colormap( &colormap_ );
  }

  void run();
  void save_eps();

  void receive( const std::string& text );
  void send( const std::string& text ) const;

protected:
  WV_Manager manager_;
  std::shared_ptr<SceneGraph> scene_;
  ClippingPlane clip_plane_;

private:
  void connect_client();
};

#ifdef AVRO_WITH_GL

class Visualizer : public Application<GLFW_Interface<OpenGL_Manager>>
{
public:
  Visualizer();

  void add_topology( const TopologyBase& topology )
  {
    // add the topology to the relevant windows
    index_t id = main_->create_scene();
    scenes_.push_back( &main_->scene(id) );
    index_t prim_id = main_->scene(id).add_primitive(topology); // create a new root in the scene graph
    manager_.select_shader( main_->scene(id).primitive(prim_id) , "wv" );
    main_->scene(id).set_colormap(&colormap_);
    for (index_t j=0;j<main_->scene(id).primitive(prim_id).nb_children();j++)
      manager_.select_shader( main_->scene(id).primitive(prim_id).child(j) , "wv" );
  }

  void remove( index_t scene , index_t root )
  {
    scenes_[scene]->remove(root);
    if (scenes_[scene]->nb_primitives()==0)
    {
      scenes_.erase( scenes_.begin() + root );
    }
  }

  OpenGL_Manager& manager() { return manager_; }

  GLFW_Window& main_window() { return *main_.get(); }

  std::shared_ptr<GLFW_Window> main_;
  //GLFW_Window side_;
private:
  std::shared_ptr<Widget> toolbar_;
};

#else

class Visualizer
{
public:
  Visualizer() {}

  void add_topology( const TopologyBase& topology )
  {}

  void run()
  {
    printf("graphics not supported");
  }
};

#endif

class WebVisualizer : public Application<Web_Interface>
{
public:

  WebVisualizer();

  void add_topology( const TopologyBase& topology )
  {
    // create a new root in the scene graph
    scene_->add_primitive(topology);
  }
private:

};

} // graphics

} // avro

#endif
