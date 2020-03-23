#ifndef avro_LIB_GRAPHICS_WINDOW_H_
#define avro_LIB_GRAPHICS_WINDOW_H_

//#include "graphics/application.h"
#include "graphics/controls.h"
#include "graphics/gl.h"
#include "graphics/math.h"
#include "graphics/user_interface.h"
#include "graphics/scene.h"

#include <memory>
#include <string>

namespace avro
{

namespace graphics
{

class GLFW_Window
{
public:
  GLFW_Window( GraphicsManager& manager , int width , int height , const std::string& title );

  void
  mouse_button_callback(int button,int action,int mods)
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
  mouse_move_callback(double xpos, double ypos)
  {
    controls_.mouse_move((int)xpos,(int)ypos);
  }

  void
  mouse_scroll_callback(double xpos, double ypos)
  {
    controls_.mouse_wheel(xpos,ypos);
  }

  void key_callback(int key, int scancode, int action, int mods)
  {
    controls_.key_down(key);
  }

  void make_current()
  {
    glfwMakeContextCurrent(window_);
  }

  void update_view();

  bool should_close()
  {
    return glfwWindowShouldClose(window_) || (glfwGetKey(window_, GLFW_KEY_ESCAPE ) == GLFW_PRESS);
  }

  void save_eps( const std::string& filename );

  void create_interface()
  {
    interface_ = std::make_shared<Interface>(*this);
  }

  void setup()
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

  void begin_draw()
  {
    make_current();

    glfwPollEvents();

    glClearColor (1.0, 1.0, 1.0, 0.0); // white
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (interface_!=nullptr)
      interface_->begin_draw();
  }

  void draw_interface() const;

  void end_draw()
  {
    if (interface_!=nullptr)
      interface_->end_draw();
    glfwSwapBuffers(window_);

    updated_ = controls_.update();
  }

  index_t nb_scene() const { return scene_.size(); }
  SceneGraph& scene( index_t k ) { return scene_[k]; }
  const SceneGraph& scene( index_t k ) const { return scene_[k]; }

  index_t create_scene()
  {
    index_t id = scene_.size();
    scene_.push_back(SceneGraph());
    return id;
  }

  Interface& interface() { return *interface_.get(); }
  const Interface& interface() const { return *interface_.get(); }

  Controls& controls() { return controls_; }

  GLFWwindow* window() { return window_; }
  const GLFWwindow* window() const { return window_; }

  int width() const { return width_; }
  int height() const { return height_; }

private:
  std::string title_;
  GLFWwindow* window_;
  GraphicsManager& manager_;

  int width_;
  int height_;
  float fov_ = 45.0f;

  vec3 position_;
  float angles_[2];

  std::vector<SceneGraph> scene_;
  Controls controls_;

  std::shared_ptr<Interface> interface_;

  bool updated_;
};

} // graphics

} // avro

#endif
