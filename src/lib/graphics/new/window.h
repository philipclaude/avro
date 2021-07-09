#ifndef AVRO_LIB_GRAPHICS_NEW_WINDOW_H_
#define AVRO_LIB_GRAPHICS_NEW_WINDOW_H_

#include "common/error.h"

#include "graphics/gl.h"

namespace avro
{

namespace graphics
{

class Window {
public:
  Window( index_t width , index_t height ) :
    width_(width),
    height_(height)
  {}

  void init() {

    avro_assert_msg( glfwInit() , "problem initializing OpenGL!" );
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 6);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

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

    //glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
  }

  GLFWwindow* window() { return window_; }

private:
  index_t width_;
  index_t height_;
  GLFWwindow* window_;

};

} // graphics

} // avro

#endif
