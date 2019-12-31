#include "graphics/controls.h"

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace avro
{

namespace graphics
{

#define   FloatInfinity std::numeric_limits<float>::infinity()
#define   SQRT1_2  0.7071067811865476

Camera::Camera(const glm::vec3& e) :
  eye(e), up(0.0f,1.0f,0.0f), viewMatrix(1.0f)
{}

void
Camera::lookAt(const glm::vec3& target)
{
  viewMatrix = glm::lookAt(eye,target,up);
}

Trackball::Trackball(Camera* cam,glm::vec4 screenSize) :
  m_rotateSpeed(1.0f),
  m_zoomSpeed(1.2f),
  m_panSpeed(0.3f),
  m_dynamicDampingFactor(0.2f),
  m_minDistance(0.0f),
  m_maxDistance(FloatInfinity),
  m_enabled(true),
  m_noRotate(false),
  m_noZoom(false),
  m_noPan(false),
  m_noRoll(false),
  m_staticMoving(false),
  camera_(cam),
  m_screen(screenSize),
  m_target(0.0f),
  m_lastPos(0.0f),
  m_eye(0.0f),
  m_rotStart(0.0f),
  m_rotEnd(0.0f),
  m_zoomStart(0.0f),
  m_zoomEnd(0.0f),
  m_panStart(0.0f),
  m_panEnd(0.0f),
  state_(TCB_STATE::NONE),
  m_prevState(TCB_STATE::NONE)
{}

void
Trackball::reset(glm::vec4 screenSize)
{
  m_screen = screenSize;
  m_enabled = true;
  m_rotateSpeed = 1.0f;
  m_zoomSpeed=1.2f;
  m_panSpeed=0.3f;
  m_noRotate=false;
  m_noPan=false;
  m_noZoom=false;
  m_noRoll=false;
  m_staticMoving=false;
  m_dynamicDampingFactor=0.2f;
  m_minDistance=0.0f;
  m_maxDistance=FloatInfinity;
  //m_target=0.0f;
  //m_lastPos=0.0f;
  state_=TCB_STATE::NONE;
  m_prevState=TCB_STATE::NONE;
  //m_eye=0.0f;
  //m_rotStart=0.0f;
  //m_rotEnd=0.0f;
  //m_zoomStart=0.0f;
  //m_zoomEnd=0.0f;
  //m_panStart=0.0f;
  //m_panEnd=0.0f;
}

void
Trackball::update()
{
  m_eye = camera_->eye - m_target;

  if (!m_noRotate)
  {
    RotateCamera();
  }

  if (!m_noZoom)
  {
    ZoomCamera();
  }

  if (!m_noPan)
  {
    PanCamera();
  }

  // object.position =  target + _eye;
  camera_->eye = (m_target + m_eye);

  CheckDistances();

  camera_->lookAt(m_target);

  if (glm::length2(m_lastPos - camera_->eye) > 0.0f)
  {
    m_lastPos = camera_->eye;
  }
}

void
Trackball::RotateCamera()
{
  float angle = (float)acos(glm::dot(m_rotStart, m_rotEnd) / glm::length(m_rotStart) / glm::length(m_rotEnd));

  if (!std::isnan(angle) && angle != 0.0f)
  {
    glm::vec3 axis = glm::normalize(glm::cross( m_rotStart, m_rotEnd ));

    if (glm::isnan(axis.x) || glm::isnan(axis.y) || glm::isnan(axis.z))
      return;

    glm::quat quaternion;

    angle *= m_rotateSpeed;

    quaternion = glm::angleAxis(-angle,axis);

    m_eye = glm::rotate( quaternion ,  m_eye);

    camera_->up =  glm::rotate(quaternion , camera_->up);

    m_rotEnd = glm::rotate( quaternion , m_rotEnd);

    if (m_staticMoving)
    {
      m_rotStart = m_rotEnd;
    }
    else
    {
      quaternion = glm::angleAxis( angle * (m_dynamicDampingFactor - 1.0f),axis);
      m_rotStart = glm::rotate(quaternion,m_rotStart);
    }
  }
}

void
Trackball::ZoomCamera()
{
  float factor = 1.0f + (float)(m_zoomEnd.y - m_zoomStart.y) * m_zoomSpeed;

  if (factor != 1.0f && factor > 0.0f)
  {
    m_eye = m_eye * (float)factor;

    if (m_staticMoving)
    {
      m_zoomStart = m_zoomEnd;
    }
    else
    {
      m_zoomStart.y += (float)(m_zoomEnd.y - m_zoomStart.y) * m_dynamicDampingFactor;
    }
  }
}

void Trackball::PanCamera()
{
  glm::vec2 mouseChange = m_panEnd - m_panStart;

  if (glm::length(mouseChange) != 0.0f)
  {
    mouseChange *= glm::length( m_eye) * m_panSpeed;

    glm::vec3 pan =glm::normalize(glm::cross(m_eye, camera_->up));

    pan *= mouseChange.x;

    glm::vec3 camUpClone = glm::normalize(camera_->up);

    camUpClone *=  mouseChange.y;
    pan += camUpClone;

    camera_->eye += pan;

    m_target +=  pan;
    if (m_staticMoving)
    {
      m_panStart = m_panEnd;
    }
    else
    {
      m_panStart += (m_panEnd - m_panStart) * m_dynamicDampingFactor;
    }
  }
}

void
Trackball::CheckDistances()
{
  if (!m_noZoom || !m_noPan)
  {
    if (glm::length2(camera_->eye)	 > m_maxDistance * m_maxDistance)
    {
      camera_->eye = glm::normalize(camera_->eye) * m_maxDistance;
    }

    if (glm::length2(m_eye) < m_minDistance * m_minDistance)
    {
      m_eye = glm::normalize(m_eye) * m_minDistance;
      camera_->eye = m_target + m_eye;
    }
  }
}

glm::vec3
Trackball::GetMouseProjectionOnBall(int clientX, int clientY)
{
  glm::vec3 mouseOnBall = glm::vec3
  (
    ((float)clientX - (float)m_screen.z * 0.5f) / (float)(m_screen.z * 0.5f),
    ((float)clientY - (float)m_screen.w * 0.5f) / (float)(m_screen.w * 0.5f),
    0.0f
  );

  double length = (double) glm::length(mouseOnBall);

  if (m_noRoll)
  {
    if (length < SQRT1_2)
    {
      mouseOnBall.z = (float)sqrt(1.0 - length * length);
    }
    else
    {
      mouseOnBall.z = (float)(0.5 / length);
    }
  }
  else if (length > 1.0)
  {
    mouseOnBall = glm::normalize(mouseOnBall);
  }
  else
  {
    mouseOnBall.z = (float)sqrt(1.0 - length * length);
  }

  m_eye = m_target - camera_->eye;

  glm::vec3 upClone = camera_->up;

  upClone = glm::normalize(upClone);
  glm::vec3 projection = upClone * mouseOnBall.y;

  glm::vec3 cross = glm::normalize( glm::cross(camera_->up,m_eye) );

  cross      *= mouseOnBall.x;
  projection += cross;

  glm::vec3 eyeClone = glm::normalize(m_eye);

  projection += eyeClone * mouseOnBall.z;

  return projection;
}

void
Trackball::MouseDown(int button, int action, int mods,int xpos,int ypos)
{
  if (!m_enabled) return;

  if (state_ == TCB_STATE::NONE)
  {
    if (button == GLFW_MOUSE_BUTTON_RIGHT)
    {
      state_ = TCB_STATE::PAN;
    }
    else
    {
      state_ = TCB_STATE::ROTATE;
    }
  }

  if (state_ == TCB_STATE::ROTATE && !m_noRotate)
  {
    m_rotStart = GetMouseProjectionOnBall(xpos, ypos);
    m_rotEnd   = m_rotStart;

  }
  else if (state_ == TCB_STATE::ZOOM && !m_noZoom)
  {
    m_zoomStart = GetMouseOnScreen(xpos, ypos);
    m_zoomEnd   = m_zoomStart;
  }
  else if (state_ == TCB_STATE::PAN && !m_noPan)
  {
    m_panStart = GetMouseOnScreen(xpos, ypos);
    m_panEnd   = m_panStart;
  }
}

void
Trackball::KeyDown(int key)
{
  if (!m_enabled) return;

  m_prevState = state_;

  if (state_ != TCB_STATE::NONE)
  {
    return;
  }
  else if ((key == GLFW_KEY_LEFT_ALT || key == GLFW_KEY_RIGHT_CONTROL) && !m_noRotate)
  {
    state_ = TCB_STATE::ROTATE;
  }
  else if (key == GLFW_KEY_Z && !m_noZoom)
  {
    state_ = TCB_STATE::ZOOM;
  }
  else if ( (key == GLFW_KEY_LEFT_SHIFT || key == GLFW_KEY_RIGHT_SHIFT) && !m_noPan)
  {
    state_ = TCB_STATE::PAN;
  }
  else if ( key == GLFW_KEY_E )
  {
    controls_.edges_visible = !controls_.edges_visible;
  }
  else if ( key == GLFW_KEY_F )
  {
    controls_.faces_visible = !controls_.faces_visible;
  }
  else if ( key == GLFW_KEY_P )
  {
    controls_.points_visible = !controls_.points_visible;
  }
  else if ( key == GLFW_KEY_R )
  {
    camera_->eye = glm::vec3(0.0f,0.0f,7.0f);
    camera_->up  = glm::vec3(0.0f,1.0f,0.0f);
    reset(glm::vec4(0.0f,0.0f,(float)1024,(float)640));
  }
}

void
Trackball::MouseWheel(double xoffset ,double yoffset)
{
  if (!m_enabled) return;

  float delta = 0.0f;

  if (yoffset != 0.0)
  {
    delta = (float) yoffset/3.0f;
  }

  m_zoomStart.y += delta*0.05f;
}

void
Trackball::MouseMove(int xpos,int ypos)
{
  if (!m_enabled) return;

  if (state_ == TCB_STATE::ROTATE && !m_noRotate)
  {
    m_rotEnd = GetMouseProjectionOnBall(xpos,ypos);
  }
  else if (state_ == TCB_STATE::ZOOM && !m_noZoom)
  {
    m_zoomEnd = GetMouseOnScreen(xpos,ypos);
  }
  else if (state_ == TCB_STATE::PAN && !m_noPan)
  {
    m_panEnd = GetMouseOnScreen(xpos,ypos);
  }
}

} // graphics

} // avro
