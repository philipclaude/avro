#ifndef AVRO_LIB_GRAPHICS_NEW_WINDOW_H_
#define AVRO_LIB_GRAPHICS_NEW_WINDOW_H_

#include "common/error.h"

#include "graphics/gl.h"
#include "graphics/new/camera.h"
#include "graphics/new/managers.h"
#include "graphics/new/plot.h"

namespace avro
{

namespace graphics
{

class Plot;

class Trackball {
public:

  Trackball() :
    dragging_(false),
    rotating_(false)
  {}

  bool dragging() const { return dragging_; }
  bool rotating() const { return rotating_; }

  void set_dragging( bool x ) { dragging_ = x; }
  void set_rotating( bool x ) { rotating_ = x; }

  void set_current_position( double x , double y ) {
    current_x_ = x;
    current_y_ = y;
  }

  void get_current_position( double& x , double& y ) {
    x = current_x_;
    y = current_y_;
  }

  mat4 get_rotation_matrix( double X , double Y ) {

    real_t X2 = X*X, Y2 = Y*Y;
    real_t q = 1 + X2 + Y2;
    real_t s = 1 - X2 - Y2;
    real_t r2 = 1/(q*q), s2 = s*s;
    real_t A = (s2 + 4*(Y2 - X2))*r2;
    real_t B = -8*X*Y*r2;
    real_t C = 4*s*X*r2;
    real_t D = (s2 + 4*(X2 - Y2))*r2;
    real_t E = 4*s*Y*r2;
    real_t F = (s2 - 4*(X2 + Y2))*r2;

    mat4 R; // initializes to zero
    R(0,0) =  A; R(1,0) =  B; R(2,0) = C;
    R(0,1) =  B; R(1,1) =  D; R(2,1) = E;
    R(0,2) = -C; R(1,2) = -E; R(2,2) = F;
    R(3,3) =  1;
    return R;
  }
  mat4 get_translation_matrix( double x , double y ) {
    mat4 T;
    T = glm::identity();

    // compute the transformation in screen space
    vec3 t = {float(x),float(y),0.0f};
    T = glm::translate( T , t );

    return T;
  }

private:
	vec3 center_;
	mat4 target_;

  bool dragging_;
  bool rotating_;

  double current_x_;
  double current_y_;
};

class Window {
public:
  Window( index_t width , index_t height ) :
    width_(width),
    height_(height),
    camera_(M_PI/4.0,width,height),
    picked_(-1)
  {}

  ~Window();

  void init();
  void compute_view();

  GLFWwindow* window() { return window_; }
  OpenGL4_Manager& manager() { return manager_; }

  void mouse_button_callback(int button, int action, int mods);
  void mouse_move_callback(double x, double y);
  void mouse_scroll_callback(double x, double y);
  void key_callback( int key , int scancode , int action , int mods );
  void resize( int width , int height);

  void add_plot( Plot* plot ) {
    plot_.push_back(plot);

    for (index_t k = 0; k < plot->nb_vao(); k++)
      manager_.write(plot->vao(k));
  }

  void draw();


private:
  index_t width_;
  index_t height_;
  GLFWwindow* window_;

  OpenGL4_Manager manager_;

	std::vector<Plot*> plot_;

  Trackball trackball_;
  Camera camera_;
  int picked_;
};

} // graphics

} // avro

#endif
