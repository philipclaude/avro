#include "common/error.h"

#include "graphics/controls.h"
#include "graphics/gl.h"
#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/primitive.h"
#include "graphics/window.h"

namespace avro
{

namespace graphics
{

static Trackball* trackball_instance = nullptr;

static void
mouse_button_callback(GLFWwindow* win ,int button,int action,int mods)
{
  if (action == GLFW_PRESS)
  {
    double xpos,ypos;
    glfwGetCursorPos(win,&xpos,&ypos);
    trackball_instance->MouseDown(button,action,mods,(int)xpos,(int)ypos);
  }
  if (action == GLFW_RELEASE)
  {
    trackball_instance->MouseUp();
  }
}

static void
mouse_move_callback(GLFWwindow* win, double xpos, double ypos)
{
  trackball_instance->MouseMove((int)xpos,(int)ypos);
}

static void
mouse_scroll_callback(GLFWwindow* win, double xpos, double ypos)
{
  trackball_instance->MouseWheel(xpos,ypos);
}

static void
keyboard_callback(GLFWwindow* win,int key,int scancode,int action,int mods)
{
  if (action == GLFW_PRESS)
  {
    trackball_instance->KeyDown(key);
  }

  if (action == GLFW_RELEASE)
  {
    trackball_instance->KeyUp();
  }
}

Window::Window( const std::string& title , Plotter* plotter) :
  title_(title),
  plotter_(plotter),
  camera_(glm::vec3(0.0f,0.0f,7.0f)),
  trackball_(&camera_,glm::vec4(0.0f,0.0f,(float)width_,(float)height_)),
  interface_(NULL)
{
  window_ = glfwCreateWindow( width_ , height_ , title_.c_str() , NULL, NULL);
  if (!window_)
  {
    glfwTerminate();
    avro_assert_not_reached;
  }
  glfwMakeContextCurrent(window_);

  trackball_instance = &trackball_;

  // initialize the trackball callbacks
  glfwSetCursorPosCallback(window_,&mouse_move_callback);
  glfwSetMouseButtonCallback(window_,&mouse_button_callback);
  glfwSetScrollCallback(window_,&mouse_scroll_callback);
  glfwSetKeyCallback(window_,&keyboard_callback);
}

void
Window::reset()
{
  camera_.eye = glm::vec3(0.0f,0.0f,7.0f);
  camera_.up  = glm::vec3(0.0f,1.0f,0.0f);
  trackball_.reset(glm::vec4(0.0f,0.0f,(float)width_,(float)height_));
}

Window::~Window()
{
  glfwDestroyWindow(window_);
}

void
Window::setMatrices()
{
  modelMatrix_ = mat4(1.0);
  projMatrix_ = glm::perspective(glm::radians(fov_), float(width_)/float(height_) , 1.0f, 100.0f);

  trackball_.update();

  // probably not the right place to put this for performance
  for (index_t k=0;k<plot_.size();k++)
  {
    plot_[k]->set_visibility( trackball_.controls() );
  }

  // compute the matrices that need to be passed to the shaders
  viewMatrix_ = camera_.viewMatrix;
  mvp_ = projMatrix_ * viewMatrix_ * modelMatrix_;
  normalMatrix_ = glm::transpose(glm::inverse(glm::mat3( modelMatrix_*viewMatrix_)));
}

void
Window::draw()
{
  glClearColor (1.0, 1.0, 1.0, 0.0); // white
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  setMatrices();

  // this is where the drawing request occurs!
  for (index_t k=0;k<plot_.size();k++)
  {
    plot_[k]->draw();
  }
}

void
Window::write()
{
  // this is where the writing to opengl buffers happens!
  for (index_t k=0;k<plot_.size();k++)
    plot_[k]->write();
}

void
Window::save( const std::string& filename )
{
  // first set the transform feedback on for all plots
  for (index_t k=0;k<plot_.size();k++)
    plot_[k]->set_transform_feedback(true);

  draw();

  // reset transform feedback
  for (index_t k=0;k<plot_.size();k++)
    plot_[k]->set_transform_feedback(false);
}

void
Window::attach( Plot_ptr plot )
{
  // temporary until more generic plots are written
  plot->setWindow(this);
  plot_.push_back( plot );
}

void
Window::set_interface( Interface* interface )
{
  interface_ = interface;
}

void
Window::run()
{
  // ensure we can capture the escape key being pressed below
  glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_TRUE);

  // hide the mouse and enable unlimited mouvement
  //glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

  glfwPollEvents();
  glfwSetCursorPos(window_, width_/2, height_/2);

  angles_[0] = 0;
  angles_[1] = 0;

  glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

  // Enable depth test
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  // Cull triangles which normal is not towards the camera
  glDisable(GL_CULL_FACE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  write();

  bool show_demo_window = true;

  while (!glfwWindowShouldClose(window_) && glfwGetKey(window_, GLFW_KEY_ESCAPE ) != GLFW_PRESS)
  {
    glfwPollEvents();

    if (interface_)
      interface_->show();

    draw();

    if (interface_)
      interface_->render();

    glfwSwapBuffers(window_);
  }

}

} // graphics

} // avro
