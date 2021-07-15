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

  Trackball() {}


private:
	vec3 center_;
	mat4 target_;
};

class Window {
public:
  Window( index_t width , index_t height ) :
    width_(width),
    height_(height),
    camera_(M_PI/4.0,width,height)
  {}

  void init() {

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
  }

  GLFWwindow* window() { return window_; }
  OpenGL4_Manager& manager() { return manager_; }

  void mouse_down();
  void mouse_up();
  void mouse_move() {

    // check if dragging

    // check if rotating/translating

    // compute the translation/rotation matrix from the trackball
    mat4 T;

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

  }
  void key_down();
  void key_up();

  void add_plot( Plot* plot );


private:
  index_t width_;
  index_t height_;
  GLFWwindow* window_;

  OpenGL4_Manager manager_;

	std::vector<Plot*> plot_;

	int picked_;

  Trackball trackball_;
  Camera camera_;

};

} // graphics

} // avro

#endif
