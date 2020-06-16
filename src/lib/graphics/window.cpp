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

#include "graphics/colormap.h"
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
  clip_controls_(fov_,width_,height_,0.1f,100.0f),
  modify_clip_plane_(false),
  show_clip_plane_(false),
  show_axes_(true),
  center_axes_(false),
  clip_plane_(3,2.5),
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

    if (!modify_clip_plane_)
      controls_.mouse_down(button,action,mods,(int)xpos,(int)ypos);
    clip_controls_.mouse_down(button,action,mods,(int)xpos,(int)ypos);
  }

  if (action == GLFW_RELEASE)
  {
    if (!modify_clip_plane_)
      controls_.mouse_up();
    clip_controls_.mouse_up();
  }
}

void
GLFW_Window::mouse_move_callback(double xpos, double ypos)
{
  if (!modify_clip_plane_)
    controls_.mouse_move((int)xpos,(int)ypos);
  clip_controls_.mouse_move((int)xpos,(int)ypos);
}

void
GLFW_Window::mouse_scroll_callback(double xpos, double ypos)
{
  //if (!modify_clip_plane_)
  controls_.mouse_wheel(xpos,ypos);
  clip_controls_.mouse_wheel(xpos,ypos);
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
GLFW_Window::write_axes()
{
  real_t len = 0.2;
  real_t x0[3] = { -0.9 , -0.9 , 0.0 };

  // option to center the axes at the true location of the scene
  if (center_axes_)
    x0[0] = x0[1] = x0[2] = 0.0;

  std::vector<real_t> points(48);
  points[ 0] = 0.0 ; points[ 1] = 0.0 ; points[ 2] = 0.0;
  points[ 3] = len ; points[ 4] = 0.0 ; points[ 5] = 0.0;
  points[ 6] = 0.0 ; points[ 7] = len ; points[ 8] = 0.0;
  points[ 9] = 0.0 ; points[10] = 0.0 ; points[11] = len;

  real_t f = 0.05;
  real_t s = 0.025;
  real_t pi = 3.141592654;

  // draw an X
  points[12] = len + f - s*std::cos(pi/3); points[13] = - s*std::sin(pi/3); points[14] = 0.0;
  points[15] = len + f + s*std::cos(pi/3); points[16] =   s*std::sin(pi/3); points[17] = 0.0;
  points[18] = len + f - s*std::cos(pi/3); points[19] =   s*std::sin(pi/3); points[20] = 0.0;
  points[21] = len + f + s*std::cos(pi/3); points[22] = - s*std::sin(pi/3); points[23] = 0.0;

  // draw a Y
  points[24] = 0.0; points[25] = len + f    ; points[26] = 0.0;
  points[27] = 0.0; points[28] = len + f + s; points[29] = 0.0;
  points[30] = -s*std::cos(pi/3); points[31] = len + f + s + s*std::sin(pi/3); points[32] = 0.0;
  points[33] =  s*std::cos(pi/3); points[34] = len + f + s + s*std::sin(pi/3); points[35] = 0.0;

  // draw a Z
  points[36] = 0.0; points[37] = -s/2; points[38] = len + f;
  points[39] = 0.0; points[40] = -s/2; points[41] = len + f + s;
  points[42] = 0.0; points[43] =  s/2; points[44] = len + f;
  points[45] = 0.0; points[46] =  s/2; points[47] = len + f + s;

  std::vector<index_t> edges(22);

  // axes edges
  edges[0] = 0; edges[1] = 1;
  edges[2] = 0; edges[3] = 2;
  edges[4] = 0; edges[5] = 3;

  // X edges
  edges[6] = 4; edges[7] = 5;
  edges[8] = 6; edges[9] = 7;

  // Y edges
  edges[10] = 8; edges[11] = 9;
  edges[12] = 9; edges[13] = 10;
  edges[14] = 9; edges[15] = 11;

  // Z edges
  edges[16] = 12; edges[17] = 13;
  edges[18] = 13; edges[19] = 14;
  edges[20] = 14; edges[21] = 15;

  for (index_t k=0;k<points.size()/3;k++)
  for (coord_t d=0;d<3;d++)
    points[3*k+d] += x0[d];

  manager_.write( "axes" , 1 , points , edges , {} , {} );
  manager_.select_shader( "axes" , "wv" );
}

void
GLFW_Window::draw_axes()
{
  if (!show_axes_) return;

  DrawingParameters params;
  params.transparency = 1.0;
  params.mvp = controls_.model_view_projection();
  params.lighting = 1.0;
  manager_.draw( "axes" , 1 , params );
}

void
GLFW_Window::draw_plane(const real_t* focus)
{
  if (!show_clip_plane_) return;

  DrawingParameters params;
  params.transparency = 0.25;
  params.mvp = clip_controls_.model_view_projection();
  params.lighting = 1.0;
  clip_plane_.plot( manager_,focus,controls_.model_view_projection() );
}

void
GLFW_Window::flip_clipping_normal()
{
  clip_plane_.flip_normal();
}

void
GLFW_Window::draw_colorbar(const Colormap& colormap, const real_t* ulim)
{
  int ncol = 10;
  std::vector<float> values;
  colormap.generate(ncol,values);

  index_t ntri = 2*(ncol - 1);
  index_t npts = 3*ntri;

  // generate the points and triangles with the colors
  std::vector<real_t> points(3*npts);
  std::vector<index_t> triangles(3*ntri);
  std::vector<real_t> colors(3*npts);
  std::vector<index_t> edges(8);

  real_t x = 2;
  real_t y = 0.;
  real_t z = 0.8;
  real_t w = 0.2;
  real_t h = 5*w;
  real_t dh = h/real_t(ncol-1);
  index_t i0,i1,i2,i3;

  index_t pt = 0;
  for (int i=0;i<ncol-1;i++)
  {
    // first point
    points[3*pt  ] = x;
    points[3*pt+1] = y +i*dh;
    points[3*pt+2] = z;
    for (index_t j=0;j<3;j++)
      colors[3*pt+j] = values[3*i+j];
    i0 = pt;
    pt++;

    // second point
    points[3*pt  ] = x +w;
    points[3*pt+1] = y +i*dh;
    points[3*pt+2] = z;
    for (index_t j=0;j<3;j++)
      colors[3*pt+j] = values[3*i+j];
    i1 = pt;
    pt++;

    // third point
    points[3*pt  ] = x;
    points[3*pt+1] = y +(i+1)*dh;
    points[3*pt+2] = z;
    for (index_t j=0;j<3;j++)
      colors[3*pt+j] = values[3*(i+1)+j];
    i2 = pt;
    pt++;

    // fourth point
    points[3*pt  ] = x +w;
    points[3*pt+1] = y +(i+1)*dh;
    points[3*pt+2] = z;
    for (index_t j=0;j<3;j++)
      colors[3*pt+j] = values[3*(i+1)+j];
    i3 = pt;
    pt++;

    triangles[6*i  ] = i0; triangles[6*i+1] = i1; triangles[6*i+2] = i2;
    triangles[6*i+3] = i1; triangles[6*i+4] = i3; triangles[6*i+5] = i2;
  }

  for (index_t k=0;k<colors.size();k++)
    colors[k] *= 256.;

  i0 = 0;
  i1 = 1;
  i2 = pt -1;
  i3 = i2 -1;

  edges[0] = i0; edges[1] = i1;
  edges[2] = i1; edges[3] = i2;
  edges[4] = i2; edges[5] = i3;
  edges[6] = i3; edges[7] = i0;

  // write to the graphics manager
  manager_.write( "colorbar" , 2 , points , edges , triangles , colors );
  manager_.select_shader( "colorbar" , "wv" );

  // draw the colorbar
  DrawingParameters params;
  params.transparency = 1.0;
  params.mvp = controls_.perspective();
  params.lighting = -1.0; // no lighting
  manager_.draw("colorbar",2,params);
}

void
GLFW_Window::end_draw()
{
  if (interface_!=nullptr)
    interface_->end_draw();
  glfwSwapBuffers(window_);

  bool clip_updated = clip_controls_.update();
  if (modify_clip_plane_ && clip_updated)
  {
    // append the rigid-body transformation set by the clip-controls to the clipping plane
    clip_plane_.append_transformation( clip_controls_.rotation() , clip_controls_.translation() );
    clip_plane_.update();
  }
  updated_ = controls_.update();
}

void
GLFW_Window::clip()
{
  for (index_t k=0;k<scene_.size();k++)
    scene_[k]->write(manager_,&clip_plane_);
}

bool&
GLFW_Window::modify_clipping_plane()
{
  return modify_clip_plane_;
}

void
GLFW_Window::reset_clip()
{
  clip_controls_.set_ui_matrix( controls_.ui_matrix() );
  clip_controls_.set_mv_matrix( controls_.model_view() );
  clip_plane_.initialize();
}

void
GLFW_Window::update_view()
{
  controls_.calculate_view();
  clip_controls_.calculate_view();

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
