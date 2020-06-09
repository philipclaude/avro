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
#include "common/tools.h"

#include "graphics/controls.h"
#include "graphics/gl.h"
#include "graphics/primitive.h"
#include "graphics/scene.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "graphics/user_interface.h"

#include "library/eps.h"

#include <fstream>

namespace avro
{

namespace graphics
{

bool
load_ppm( const std::string& filename , GLFWimage& image )
{
  unsigned int height;
  unsigned int width;
  unsigned int max_col_value;
  unsigned int size;

  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);
  if (!file.is_open())
  {
    printf("could not open file %s\n",filename.c_str());
    return false;
  }

  // make sure we are reading a P6 ppm image
  std::string line;
  std::getline(file,line);
  if (line != "P6")
  {
    printf("unrecognized file format: %s (expected P6)\n",line.c_str());
    return false;
  }

  // skip any line that start with #
  std::getline(file,line);
  while (line[0] == '#')
    std::getline(file,line);

  // get the dimensions of the image
  std::stringstream dimensions(line);
  try
  {
    dimensions >> width;
    dimensions >> height;
  }
  catch (std::exception &e)
  {
    printf("header format error: %s\n",e.what());
    return false;
  }

  // get the maximum color value (usually 255)
  std::getline(file,line);
  std::stringstream max_val(line);
  try
  {
    max_val >> max_col_value;
  }
  catch (std::exception &e)
  {
    printf("header format error: %s\n",e.what());
    return false;
  }

  // allocate the pixels
  size = width*height;
  image.pixels = new unsigned char[3*size];
  image.width  = width;
  image.height = height;

  // save the pixels
  char aux;
  unsigned char r,g,b;
  for (index_t i = 0; i < size; i++)
  {
    file.read(&aux, 1);
    r = (unsigned char) aux;
    file.read(&aux, 1);
    g = (unsigned char) aux;
    file.read(&aux, 1);
    b = (unsigned char) aux;

    image.pixels[3*i  ] = r;
    image.pixels[3*i+1] = g;
    image.pixels[3*i+2] = b;
  }

  file.close();

  return true;
}

static void
_mouse_button_callback(GLFWwindow* window ,int button,int action,int mods)
{
  static_cast<GLFW_Window*>(glfwGetWindowUserPointer(window))->mouse_button_callback(button,action,mods);
}

static void
_mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
  static_cast<GLFW_Window*>(glfwGetWindowUserPointer(window))->mouse_move_callback(xpos,ypos);
}

static void
_mouse_scroll_callback(GLFWwindow* window, double xpos, double ypos)
{
  static_cast<GLFW_Window*>(glfwGetWindowUserPointer(window))->mouse_scroll_callback(xpos,ypos);
}

static void
_keyboard_callback(GLFWwindow* window,int key,int scancode,int action,int mods)
{
  static_cast<GLFW_Window*>(glfwGetWindowUserPointer(window))->key_callback(key,scancode,action,mods);
}

GLFW_Window::GLFW_Window( GraphicsManager& manager , int width , int height , const std::string& title ) :
  title_(title),
  manager_(manager),
  width_(width),
  height_(height),
  controls_(fov_,width_,height_,0.1f,100.0f),
  interface_(nullptr),
  updated_(true),
  fps_(60)
{
  window_ = glfwCreateWindow( width_ , height_ , title_.c_str() , NULL, NULL);
  if (!window_)
  {
    // a common error is for the graphics drivers to not support
    // the request opengl core profile (e.g. using 4.1 when drivers support 3.3)
    const char* description;
    int code = glfwGetError(&description);
    UNUSED(code);

    if (description)
      printf("GLFW error: %s\n",description);

    glfwTerminate();
    avro_assert_not_reached;
  }
  glfwMakeContextCurrent(window_);

  // save the window for performing callbacks with the correct trackball
  glfwSetWindowUserPointer(window_, this);

  // initialize the trackball callbacks
  glfwSetCursorPosCallback(window_,&_mouse_move_callback);
  glfwSetMouseButtonCallback(window_,&_mouse_button_callback);
  glfwSetScrollCallback(window_,&_mouse_scroll_callback);
  glfwSetKeyCallback( window_ , &_keyboard_callback );

  // this doesn't work on mac because you have to bundle things :(
  if (load_ppm( std::string(AVRO_SOURCE_DIR)+"/doc/fig/avro.ppm", images_[0] ))
  {
    printf("--> loaded %u x %u icon\n",images_[0].width,images_[0].height);
    glfwSetWindowIcon(window_,1,images_);
    delete [] images_[0].pixels; // free up the memory used by the pixels
  }
}

void
GLFW_Window::mouse_button_callback(int button,int action,int mods)
{
  if (action == GLFW_PRESS)
  {
    double xpos,ypos;
    glfwGetCursorPos(window_,&xpos,&ypos);
    controls_.mouse_down(button,action,mods,(int)xpos,(int)ypos);
  }
  if (action == GLFW_RELEASE)
  {
    controls_.mouse_up();
  }
}

void
GLFW_Window::mouse_move_callback(double xpos, double ypos)
{
  controls_.mouse_move((int)xpos,(int)ypos);
}

void
GLFW_Window::mouse_scroll_callback(double xpos, double ypos)
{
  controls_.mouse_wheel(xpos,ypos);
}

void
GLFW_Window::key_callback(int key, int scancode, int action, int mods)
{
  controls_.key_down(key);
}

void
GLFW_Window::make_current()
{
  glfwMakeContextCurrent(window_);
}

void
GLFW_Window::create_interface()
{
  interface_ = std::make_shared<Interface>(*this,manager_.listener());
}

bool
GLFW_Window::should_close()
{
  return glfwWindowShouldClose(window_) || (glfwGetKey(window_, GLFW_KEY_ESCAPE ) == GLFW_PRESS);
}

void
GLFW_Window::setup()
{
  // ensure we can capture the escape key being pressed below
  glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_TRUE);

  // hide the mouse and enable unlimited mouvement
  //glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

  glfwPollEvents();
  glfwSetCursorPos(window_, width_/2, height_/2);

  glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

  // enable depth test
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  // option to cull triangles which normal is not towards the camera
  glDisable(GL_CULL_FACE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1.5);
}

void
GLFW_Window::begin_draw()
{
  make_current();

  glfwPollEvents();

  glClearColor (1.0, 1.0, 1.0, 0.0); // white
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (interface_!=nullptr)
    interface_->begin_draw();
}

void
GLFW_Window::end_draw()
{
  if (interface_!=nullptr)
    interface_->end_draw();
  glfwSwapBuffers(window_);

  updated_ = controls_.update();
}

void
GLFW_Window::update_view()
{
  controls_.calculate_view();
  for (index_t k=0;k<scene_.size();k++)
  {
    scene_[k]->update_matrices(controls_);
    scene_[k]->set_update(true);
  }
}

void
GLFW_Window::save_eps( const std::string& filename )
{
  TransformFeedbackResult feedback;
  for (index_t k=0;k<scene_.size();k++)
  {
    manager_.draw(*scene_[k].get(),&feedback);
  }

  library::epsFile eps;
  int viewport[4] = {0,0,width_,height_};
  eps.set_viewport(viewport);
  eps.add_triangles( feedback.triangle_points() , feedback.triangle_colors() );
  eps.add_edges( feedback.edge_points() , feedback.edge_colors() );
  eps.write( filename );
}

index_t
GLFW_Window::create_scene()
{
  index_t id = scene_.size();
  scene_.push_back(std::make_shared<SceneGraph>());
  return id;
}

} // graphics

} // avro
