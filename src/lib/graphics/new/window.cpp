#include "graphics/new/window.h"
#include "graphics/new/shader_library.h"

namespace avro
{

namespace graphics
{

static void
_mouse_button_callback(GLFWwindow* window ,int button,int action,int mods) {
  static_cast<Window*>(glfwGetWindowUserPointer(window))->mouse_button_callback(button,action,mods);
}

static void
_mouse_move_callback(GLFWwindow* window, double x, double y) {
  static_cast<Window*>(glfwGetWindowUserPointer(window))->mouse_move_callback(x,y);
}

static void
_mouse_scroll_callback(GLFWwindow* window, double x, double y) {
  static_cast<Window*>(glfwGetWindowUserPointer(window))->mouse_scroll_callback(x,y);
}

static void
_keyboard_callback(GLFWwindow* window,int key,int scancode,int action,int mods) {
  static_cast<Window*>(glfwGetWindowUserPointer(window))->key_callback(key,scancode,action,mods);
}

static void
_resize_callback( GLFWwindow* window, int width, int height) {
  static_cast<Window*>(glfwGetWindowUserPointer(window))->resize(width,height);
}

void
Window::init() {

  window_ = glfwCreateWindow( width_ , height_ , "avro" , NULL, NULL);
  if (!window_) {
    const char* description;
    int code = glfwGetError(&description);

    if (description)
      printf("GLFW error (%d): %s\n",code,description);

    glfwTerminate();
    avro_assert_not_reached;
  }
  glfwMakeContextCurrent(window_);
  glfwSwapInterval(0); // otherwise rendering lags a bit behind cursor movement

  // load GL functions
  gladLoadGL();
  dumpGLInfo();

  // initialize the shaders (must happen after the GL functions are loaded)
  __shaders__ = std::make_shared<Shaders>(0,3,1,2);

  // ensure we can capture the escape key being pressed
  glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_TRUE);
  glfwPollEvents();
  glfwSetCursorPos(window_, width_/2, height_/2);
  glfwSetWindowSize(window_,width_,height_);
  glViewport(0, 0, width_,height_);

  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_PROGRAM_POINT_SIZE);

  // save the window for performing callbacks with the correct trackball
  glfwSetWindowUserPointer(window_, this);

  // initialize the trackball callbacks
  glfwSetCursorPosCallback(window_,&_mouse_move_callback);
  glfwSetMouseButtonCallback(window_,&_mouse_button_callback);
  glfwSetScrollCallback(window_,&_mouse_scroll_callback);
  glfwSetKeyCallback( window_ , &_keyboard_callback );
  glfwSetWindowSizeCallback( window_ , &_resize_callback );
}

void
Window::mouse_button_callback(int button, int action, int modifier) {

  if (action == GLFW_PRESS) {
    double x,y;
    glfwGetCursorPos(window_,&x,&y);
    trackball_.set_current_position(x,y);
    trackball_.set_dragging(true);

    trackball_.set_rotating(true);
    if (modifier == GLFW_MOD_SHIFT || modifier == GLFW_MOD_CONTROL)
      trackball_.set_rotating(false);
  }
  else if (action == GLFW_RELEASE) {
    trackball_.set_dragging(false);
  }
}

void
Window::mouse_move_callback(double x, double y) {

  // check if dragging
  if (!trackball_.dragging()) return;

  real_t xm,ym;
  trackball_.get_current_position(xm,ym);

  // check if rotating/translating
  // compute the translation/rotation matrix from the trackball
  mat4 T;
  if (trackball_.rotating()) {
    real_t dx = (xm - x)/width_;
    real_t dy = (ym - y)/height_;
    T = trackball_.get_rotation_matrix(dx,dy);
  }
  else {
    real_t dx = (xm - x)/width_;
    real_t dy = (ym - y)/height_;
    T = trackball_.get_translation_matrix(dx,dy);
  }

  // apply the transformation to the plots (or a single picked plot)
  if (picked_ < 0) {
    // apply the transformation to each plot's model matrix
    for (index_t k = 0; k < plot_.size(); k++) {
      plot_[k]->transform(T,true);
      if (!trackball_.rotating()) plot_[k]->transform_center(T);
    }
  }
  else {
    // apply the transformation to the particular object
    plot_[picked_]->transform(T,true);
    if (!trackball_.rotating()) plot_[picked_]->transform_center(T);
  }

  trackball_.set_current_position(x,y);
  draw();
}

void
Window::mouse_scroll_callback(double dx, double dy) {

  // move the camera eye according to the scroll
  vec3 v;
  for (coord_t d = 0; d < 3; d++)
    v(d) = camera_.eye()(d) + 0.1*dy*(camera_.lookat()(d) - camera_.eye()(d));
  camera_.set_eye(v);
  draw();
}

void
Window::key_callback( int key , int scancode , int action , int mods ) {
  // no keys are currently implemented
}

void
Window::compute_view() {

  // get each plot's center and take the average
  vec3 center;
  center.zero();

  float d = 0.0;
  for (index_t k = 0; k < plot_.size(); k++) {
    for (coord_t d = 0; d < 3; d++)
      center(d) = center(d) + plot_[k]->center()(d);
    if (plot_[k]->length_scale() > d)
      d = plot_[k]->length_scale();
  }

  for (coord_t d = 0; d < 3; d++)
    center(d) /= plot_.size();
  camera_.set_lookat(center);

  // set the camera to be in the -z direction from the center (up is +y)
  vec3 eye = center;
  eye(2) = eye(2) - 4*d; // multiply the distance a bit so we're not too close to the plots
  camera_.set_eye(eye);
}

void
Window::resize(int width, int height) {
  width_  = width;
  height_ = height;
  camera_.compute_projection(width_,height_);
  glViewport(0, 0, width_,height_);
  draw();
}

void
Window::draw() {

  float col = 1.0;
  glClearColor(col,col,col, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  const mat4& view_matrix = camera_.view_matrix();
  const mat4& projection_matrix = camera_.projection_matrix();

  for (index_t k = 0; k < plot_.size(); k++) {

    // calculate the matrices
    const mat4& model_matrix = plot_[k]->model_matrix();

    // retrieve the current vao and draw
    VertexAttributeObject& vao = plot_[k]->active_vao();
    vao.draw(model_matrix,view_matrix,projection_matrix);
  }

  glfwSwapBuffers(window_);
}

Window::~Window() {
  glfwDestroyWindow(window_);
}

} // graphics

} // avro
