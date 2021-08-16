#include "graphics/gui.h"
#include "graphics/vao.h"
#include "graphics/window.h"
#include "graphics/shader_library.h"

#include <imgui/imgui.h>
#include <json/json.hpp>

#include <fstream>

namespace avro
{

namespace graphics
{

#if AVRO_WITH_GL

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

#endif

void
Window::init() {

  #if AVRO_WITH_GL

  #if AVRO_HEADLESS_GRAPHICS
  glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
  #endif
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
  glfwSwapInterval(1); // otherwise rendering lags a bit behind cursor movement

  // load GL functions
  gladLoadGL();
  dumpGLInfo();

  // initialize the shaders (must happen after the GL functions are loaded)
  __shaders__ = std::make_shared<Shaders>(-1,3,1,2);

  // ensure we can capture the escape key being pressed
  glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_TRUE);
  glfwPollEvents();
  glfwSetCursorPos(window_, width_/2, height_/2);
  glfwSetWindowSize(window_,width_,height_);

  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_PROGRAM_POINT_SIZE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_POLYGON_OFFSET_FILL);
  glPolygonOffset(1.0, 1.5);

  // save the window for performing callbacks with the correct trackball
  glfwSetWindowUserPointer(window_, this);

  // initialize the trackball callbacks
  glfwSetCursorPosCallback(window_,&_mouse_move_callback);
  glfwSetMouseButtonCallback(window_,&_mouse_button_callback);
  glfwSetScrollCallback(window_,&_mouse_scroll_callback);
  glfwSetKeyCallback( window_ , &_keyboard_callback );
  glfwSetWindowSizeCallback( window_ , &_resize_callback );

  picking_  = false;
  picked_   = 0;

  needs_drawing_ = true;
  draw_count_ = 0;

  select_colormap("bwr");

  #else
  avro_assert_not_reached;
  #endif

  screen_matrix_(0,0) = width_/2.0;
  screen_matrix_(1,1) = height_/2.0;
  screen_matrix_(0,3) = (width_-1)/2.0;
  screen_matrix_(1,3) = (height_-1)/2.0;
  screen_matrix_(2,2) = 1.0;
  screen_matrix_(3,3) = 1.0;
}

void
Window::mouse_button_callback(int button, int action, int modifier) {

  if (!enable_controls_) return;

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

  if (!enable_controls_) return;

  // check if dragging
  if (!trackball_.dragging() && !picking_) return;

  if (picking_) {
    // compute the world coordinates of the mouse
    mat4 ms = glm::identity();

    ms(0,0) =  (float)width_/2.;
    ms(1,1) =  (float)height_/2.;
    ms(0,3) =  (float)width_/2.;
    ms(1,3) =  (float)height_/2.;

    mat4 cs = ms * camera_.projection_matrix() * camera_.view_matrix();

    // find the closest plot
    real_t d = 1e20;
    for (index_t k = 0; k < plot_.size(); k++) {
      mat4 mcs = cs * plot_[k]->model_matrix();

      // transform the plot center to screen coordinates
      vec4 c = mcs * glm::to_vec4( plot_[k]->center() , 1.0 );

      real_t dx = (x - c(0)/c(3));
      real_t dy = (y - c(1)/c(3));
      real_t distance = dx*dx + dy*dy;
      if (distance < d) {
        d        = distance;
        picked_  = k;
      }
    }
  }

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
  if (!picking_) {
    // apply the transformation to each plot's model matrix
    for (index_t k = 0; k < plot_.size(); k++) {
      plot_[k]->transform(T,false);
      if (!trackball_.rotating()) plot_[k]->transform_center(T);
    }
  }
  else {
    // apply the transformation to the particular object
    plot_[picked_]->transform(T,true);
    if (!trackball_.rotating()) plot_[picked_]->transform_center(T);
  }

  trackball_.set_current_position(x,y);
  needs_drawing_ = true;
}

void
Window::mouse_scroll_callback(double dx, double dy) {

  if (!enable_controls_) return;

  // move the camera eye according to the scroll
  vec3 v;
  for (coord_t d = 0; d < 3; d++)
    v(d) = camera_.eye()(d) + 0.1*dy*(camera_.lookat()(d) - camera_.eye()(d));
  camera_.set_eye(v);
  needs_drawing_ = true;
}

void
Window::key_callback( int key , int scancode , int action , int mods ) {

  if (!enable_controls_) return;

  if (key == GLFW_KEY_P) {
    if (action == GLFW_PRESS) {
      picking_ = !picking_;
    }
    printf("picking is %s\n",picking_? "on" : "off");
    picked_ = -1;
  }
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
Window::compute_view( const vec3& center , float d ) {

  camera_.set_lookat(center);

  // set the camera to be in the -z direction from the center (up is +y)
  vec3 eye = center;
  eye(2) = eye(2) - d; // multiply the distance a bit so we're not too close to the plots
  camera_.set_eye(eye);
}

void
Window::resize(int width, int height) {
  width_  = width;
  height_ = height;
  camera_.compute_projection(width_,height_);
  needs_drawing_ = true;
  #if __APPLE__
  #else
  glViewport(0,0,width_,height_);
  #endif
  glfwSetWindowSize(window_,width_,height_);
}

void
Window::draw(bool swap_buffers) {

  if (draw_callback_ != nullptr) {
    draw_callback_(*this);
    return;
  }

  // only draw if we need to
  if (!needs_drawing_) return;
  draw_count_++;

  float col = 1.0;
  glClearColor(col,col,col, 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  const mat4& view_matrix = camera_.view_matrix();
  const mat4& projection_matrix = camera_.projection_matrix();

  for (index_t k = 0; k < plot_.size(); k++) {

    // retrieve the model matrix from the plot
    const mat4& model_matrix = plot_[k]->model_matrix();

    plot_[k]->clip().draw(view_matrix,projection_matrix);

    // retrieve the current vao and draw
    VertexAttributeObject& vao = plot_[k]->active_vao();
    const ClipPlane& clip = plot_[k]->clip();
    vao.set_lighting(lighting_);
    vao.draw(model_matrix,view_matrix,projection_matrix,&clip);
  }

  if (swap_buffers) {
    glfwSwapBuffers(window_);
  }
  needs_drawing_ = false;
}

void
Window::load_view( const std::string& filename ) {

  nlohmann::json jv;
  std::ifstream file_in(filename);
  file_in >> jv;

  float fov = jv["fov"];
  float width = jv["width"];
  float height = jv["height"];

  std::vector<float> eye_data = jv["eye"];
  std::vector<float> lookat_data = jv["lookat"];
  std::vector<float> model_matrix_data = jv["model_matrix"];

  vec3 eye(eye_data);
  vec3 lookat(lookat_data);
  mat4 model_matrix(model_matrix_data);

  resize(width,height);
  camera_.set_fov(fov);
  camera_.set_eye(eye);
  camera_.set_lookat(lookat);
  camera_.compute_projection(width,height);

  // set each plot to have the same model matrix
  // (this way we don't need to make sure the same scene is loaded)
  for (index_t k = 0; k < nb_plots(); k++)
    plot_[k]->set_model_matrix(model_matrix);

}

void
Window::save_view( const std::string& filename ) {

  nlohmann::json jv;

  // only grab the first model matrix
  avro_assert( nb_plots() > 0 );
  const mat4& model_matrix = plot_[0]->model_matrix();
  std::vector<float> model_matrix_data( model_matrix.data() , model_matrix.data()+16 );
  jv["model_matrix"] = model_matrix_data;

  std::vector<float> eye( camera_.eye().data() , camera_.eye().data()+3 );
  std::vector<float> lookat( camera_.lookat().data() , camera_.lookat().data()+3 );
  jv["eye"] = eye;
  jv["lookat"] = lookat;

  jv["fov"] = camera_.fov();
  jv["width"] = width_;
  jv["height"] = height_;

  std::ofstream file(filename);
  file << jv;
  file.close();
}

Window::~Window() {
  glfwDestroyWindow(window_);
}

} // graphics

} // avro
