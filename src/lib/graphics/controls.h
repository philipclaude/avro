#ifndef AVRO_LIB_GRAPHICS_CONTROLS_H_
#define AVRO_LIB_GRAPHICS_CONTROLS_H_

#include "graphics/gl.h"
#include "graphics/math.h"

namespace avro
{

namespace graphics
{

class Controls
{

public:

  Controls( float fov , int width , int height , float znear=0.1f , float zfar=100.0f );

  bool update();
  void calculate_view();

  void mouse_down(int button, int action, int mods,int xpos,int ypos);
  void mouse_up();
  void mouse_move(int x, int y);
  void key_down(int key);
  void mouse_wheel(double xoffset ,double yoffset);

  const glm::mat4& model_view_projection() const {return model_view_projection_; }
  const glm::mat4& normal() const { return normal_; }

  bool dragging;
  int modifier;

  void disable() { enabled_ = false; }
  void enable() { enabled_ = true; }

private:
  glm::mat4 perspective_;
  glm::mat4 model_view_;
  glm::mat4 model_view_projection_;
  glm::mat4 normal_;
  glm::mat4 ui_matrix_;

  int width_;
  int height_;

  int offleft_;
  int offtop_;

  int startx_;
  int starty_;

  int cursorx_;
  int cursory_;

  float scale_;

  const glm::vec3 eye_ = {0,0,7};
  const glm::vec3 center_ = {0,0,0};
  const glm::vec3 up_ = {0,1,0};

  bool enabled_;
};

} // graphics

} // avro

#endif
