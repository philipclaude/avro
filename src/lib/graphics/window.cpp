#include "common/error.h"

#include "graphics/controls.h"
#include "graphics/gl.h"
#include "graphics/primitive.h"
#include "graphics/window.h"

#include "graphics/user_interface.h"

#include "library/eps.h"

namespace avro
{

namespace graphics
{

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
  manager_(manager),
  title_(title),
  width_(width),
  height_(height),
  //camera_(glm::vec3(0.0f,0.0f,7.0f)),
  //trackball_(&camera_,glm::vec4(0.0f,0.0f,(float)width_,(float)height_),this),
  interface_(nullptr),
  controls_(45.0f,width_,height_,0.1f,100.0f),
  updated_(true)
{
  window_ = glfwCreateWindow( width_ , height_ , title_.c_str() , NULL, NULL);
  if (!window_)
  {
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
}

void
GLFW_Window::update_view()
{
  controls_.calculate_view();
  for (index_t k=0;k<scene_.size();k++)
  {
    scene_[k].update_matrices(controls_);
    scene_[k].set_update(true);
  }
}

void
GLFW_Window::save_eps( const std::string& filename )
{
  TransformFeedbackResult feedback;
  for (index_t k=0;k<scene_.size();k++)
  {
    printf("capture with transform feedback!!\n");
    manager_.draw(scene_[k],&feedback);
  }


  library::epsFile eps;
  int viewport[4] = {0,0,1024,640};
  eps.set_viewport(viewport);
  eps.add_triangles( feedback.triangle_points() , feedback.triangle_colors() );
  eps.add_edges( feedback.edge_points() , feedback.edge_colors() );
  eps.write( filename );
}

} // graphics

} // avro
