#include "graphics/new/window.h"

namespace avro
{

namespace graphics
{

static void
_mouse_button_callback(GLFWwindow* window ,int button,int action,int mods)
{
  static_cast<Window*>(glfwGetWindowUserPointer(window))->mouse_button_callback(button,action,mods);
}

static void
_mouse_move_callback(GLFWwindow* window, double x, double y)
{
  static_cast<Window*>(glfwGetWindowUserPointer(window))->mouse_move_callback(x,y);
}

static void
_mouse_scroll_callback(GLFWwindow* window, double x, double y)
{
  static_cast<Window*>(glfwGetWindowUserPointer(window))->mouse_scroll_callback(x,y);
}

static void
_keyboard_callback(GLFWwindow* window,int key,int scancode,int action,int mods)
{
  static_cast<Window*>(glfwGetWindowUserPointer(window))->key_callback(key,scancode,action,mods);
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

  // ensure we can capture the escape key being pressed
  glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_TRUE);
  glfwPollEvents();
  glfwSetCursorPos(window_, width_/2, height_/2);

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

  // check if rotating/translating
  // compute the translation/rotation matrix from the trackball
  mat4 T;
  if (trackball_.rotating()) {
    real_t xm,ym;
    trackball_.get_current_position(xm,ym);
    real_t dx = -(x - xm)/width_;
    real_t dy =  (y - ym)/height_;
    T = trackball_.get_rotation_matrix(dx,dy);
  }
  else {
    T = trackball_.get_translation_matrix(x,y);
  }

  // apply the transformation to the plots (or a single picked plot)
  if (picked_ < 0) {
    // apply the transformation to each plot's model matrix
    for (index_t k = 0; k < plot_.size(); k++) {
      plot_[k]->transform(T,false);
    }
  }
  else {
    // apply the transformation to the particular object
    plot_[picked_]->transform(T,true);
  }

  trackball_.set_current_position(x,y);
  draw();

}

void
Window::mouse_scroll_callback(double dx, double dy) {

  // move the camera or scale the model matrix?

  printf("dy = %g\n",dy);
  vec3 v;
  for (coord_t d = 0; d < 3; d++)
    v(d) = camera_.eye()(d) + 0.001*dy*(camera_.lookat()(d) - camera_.eye()(d));

  camera_.set_lookat(v);
  v.print();

  draw();
}

void
Window::key_callback( int key , int scancode , int action , int mods ) {
  // no keys are currently implemented
}

void
Window::draw() {

  float col = 1.0;
  glClearColor(col,col,col, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  const mat4& view_matrix = camera_.view_matrix();
  const mat4& perspective_matrix = camera_.projection_matrix();

  for (index_t k = 0; k < plot_.size(); k++) {

    // calculate the matrices
    const mat4& model_matrix = plot_[k]->model_matrix();

    mat4 mv  = view_matrix * model_matrix;
    mat4 mvp = perspective_matrix * mv;
    mat4 normal_matrix = glm::transpose(glm::inverse(mv));

    ShaderProgram* ts = plot_[k]->triangle_shader();
    ts->use();
    ts->setUniform("u_ModelViewProjectionMatrix",mvp);
    ts->setUniform("u_NormalMatrix",normal_matrix);
    ts->setUniform("u_ModelViewMatrix",mv);

    ShaderProgram* es = plot_[k]->edge_shader();
    es->use();
    es->setUniform("u_ModelViewProjectionMatrix",mvp);

    ShaderProgram* ps = plot_[k]->point_shader();
    ps->use();
    ps->setUniform("u_ModelViewProjectionMatrix",mvp);

    // retrieve the current vao
    VertexAttributeObject& vao = plot_[k]->active_vao();

    // draw the primitives
    vao.draw_edges(*es);
    vao.draw_triangles(*ts);
    //vao.draw_points(*ps);
  }

  glfwSwapBuffers(window_);
}

} // graphics

} // avro
