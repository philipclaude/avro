#ifndef URSA_LIB_GRAPHICS_CONTROLS_H_
#define URSA_LIB_GRAPHICS_CONTROLS_H_

#include "graphics/gl.h"
#include "graphics/math.h"

namespace ursa
{

namespace graphics
{

class Camera
{
public:
  Camera(const glm::vec3& pos);

  void lookAt(const glm::vec3& target);

  glm::vec3 eye, up;
  glm::mat4 viewMatrix;
  glm::mat4 mvp_;
  glm::mat4 mv_;
};

class Trackball
{
public:
  Trackball(Camera* cam,glm::vec4 screenSize);

  void update();
  void MouseDown(int button, int action, int mods,int xpos,int ypos);
  void MouseMove(int xpos, int ypos);
  void KeyDown(int key);
  void MouseWheel(double xoffset ,double yoffset);

  float m_rotateSpeed;
  float m_zoomSpeed;
  float m_panSpeed;
  float m_dynamicDampingFactor;
  float m_minDistance;
  float m_maxDistance;
  bool m_enabled;
  bool m_noRotate;
  bool m_noZoom;
  bool m_noPan;
  bool m_noRoll;
  bool m_staticMoving;

  inline void MouseUp()
  {
    if (!m_enabled) return;
    state_ = TCB_STATE::NONE;
  }

  inline void KeyUp()
  {
    if (!m_enabled) return;
    state_ = m_prevState;
  }

private:

  glm::vec3 GetMouseProjectionOnBall(int clientX, int clientY);

  inline glm::vec2
  GetMouseOnScreen(int clientX, int clientY)
  {
    return glm::vec2(
      (float)(clientX - m_screen.x) / m_screen.z,
      (float)(clientY - m_screen.y) / m_screen.w
      );
  }

  void RotateCamera();
  void ZoomCamera();
  void PanCamera();
  void CheckDistances();

  enum class TCB_STATE:uint8_t
  {
    NONE,
    ROTATE,
    ZOOM,
    PAN
  };

  Camera*   camera_;
  glm::vec4 m_screen;

  glm::vec3 m_target;
  glm::vec3 m_lastPos;
  glm::vec3 m_eye;
  glm::vec3 m_rotStart;
  glm::vec3 m_rotEnd;
  glm::vec2 m_zoomStart;
  glm::vec2 m_zoomEnd;
  glm::vec2 m_panStart;
  glm::vec2 m_panEnd;
  TCB_STATE state_;
  TCB_STATE m_prevState;
};

} // graphics

} // ursa

#endif
