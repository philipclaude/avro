#include "common/error.h"

#include "graphics/controls.h"
#include "graphics/gl.h"
#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"


namespace ursa
{

namespace graphics
{

static Trackball* trackball_instance = nullptr;

static void
mouse_button_callback(GLFWwindow* win ,int button,int action,int mods)
{
  if (action == GLFW_PRESS)
  {
    double xpos,ypos;
    glfwGetCursorPos(win,&xpos,&ypos);
    trackball_instance->MouseDown(button,action,mods,(int)xpos,(int)ypos);
  }
  if (action == GLFW_RELEASE)
  {
    trackball_instance->MouseUp();
  }
}

static void
mouse_move_callback(GLFWwindow* win, double xpos, double ypos)
{
  trackball_instance->MouseMove((int)xpos,(int)ypos);
}

static void
mouse_scroll_callback(GLFWwindow* win, double xpos, double ypos)
{
  trackball_instance->MouseWheel(xpos,ypos);
}

static void
keyboard_callback(GLFWwindow* win,int key,int scancode,int action,int mods)
{
  if (action == GLFW_PRESS)
  {
    trackball_instance->KeyDown(key);
  }

  if (action == GLFW_RELEASE)
  {
    trackball_instance->KeyUp();
  }
}

Window::Window( const std::string& title , Plotter* plotter) :
  title_(title),
  plotter_(plotter),
  camera_(glm::vec3(0.0f,0.0f,100.0f)),
  trackball_(&camera_,glm::vec4(0.0f,0.0f,(float)width_,(float)height_))
{
  window_ = glfwCreateWindow( width_ , height_ , title_.c_str() , NULL, NULL);
  if (!window_)
  {
    glfwTerminate();
    ursa_assert_not_reached;
  }
  glfwMakeContextCurrent(window_);

  trackball_instance = &trackball_;

  // initialize the trackball callbacks
  glfwSetCursorPosCallback(window_,&mouse_move_callback);
  glfwSetMouseButtonCallback(window_,&mouse_button_callback);
  glfwSetScrollCallback(window_,&mouse_scroll_callback);
  glfwSetKeyCallback(window_,&keyboard_callback);

}

Window::~Window()
{
  glfwDestroyWindow(window_);
}

void
Window::setMatrices()
{
  modelMatrix_ = mat4(1.0);
  // Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
  projMatrix_ = glm::perspective(glm::radians(fov_), 4.0f / 3.0f, 0.1f, 100.0f);

  #if 0

  // glfwGetTime is called only once, the first time this function is called
	static double lastTime = glfwGetTime();

	// Compute time difference between current and last frame
	double currentTime = glfwGetTime();
	float deltaTime = float(currentTime - lastTime);

	// Get mouse position
	double xpos, ypos;
	glfwGetCursorPos(window_, &xpos, &ypos);

	// Reset mouse position for next frame
	glfwSetCursorPos(window_, width_/2, height_/2);

	// Compute new orientation
	angles_[0] += mouseSpeed_ * float(width_/2  - xpos );
	angles_[1] += mouseSpeed_ * float(height_/2 - ypos );

	// Direction : Spherical coordinates to Cartesian coordinates conversion
	glm::vec3 direction(
		cos(angles_[1]) * sin(angles_[0]),
		sin(angles_[1]),
		cos(angles_[1]) * cos(angles_[0]) );

	// Right vector
	glm::vec3 right = glm::vec3(
		sin(angles_[0] - 3.14f/2.0f),
		0,
		cos(angles_[0] - 3.14f/2.0f) );

	// Up vector
	glm::vec3 up = glm::cross( right, direction );

	// Move forward
	if (glfwGetKey( window_, GLFW_KEY_UP ) == GLFW_PRESS)
  {
		position_ += direction * deltaTime * speed_;
	}

	// Move backward
	if (glfwGetKey( window_, GLFW_KEY_DOWN ) == GLFW_PRESS)
  {
		position_ -= direction * deltaTime * speed_;
	}

	// Strafe right
	if (glfwGetKey( window_, GLFW_KEY_RIGHT ) == GLFW_PRESS)
  {
		position_ += right * deltaTime * speed_;
	}

	// Strafe left
	if (glfwGetKey( window_, GLFW_KEY_LEFT ) == GLFW_PRESS)
  {
		position_ -= right * deltaTime * speed_;
	}

  //fov_ = fov_ - 5 * glfwGetMouseWheel();

	// Camera matrix
	viewMatrix_ = glm::lookAt(
								position_,           // Camera is here
								position_+direction, // and looks here : at the same position, plus "direction"
								up                  // Head is up (set to 0,-1,0 to look upside-down)
						   );

  mvp_ = projMatrix_*viewMatrix_*modelMatrix_;

  normalMatrix_ = glm::transpose(glm::inverse(glm::mat3( modelMatrix_*viewMatrix_)));

	// For the next frame, the "last time" will be "now"
	lastTime = currentTime;
  #else

  trackball_.update();

  viewMatrix_ = camera_.viewMatrix;
  mvp_ = projMatrix_ * viewMatrix_ * modelMatrix_;
  normalMatrix_ = glm::transpose(glm::inverse(glm::mat3( modelMatrix_*viewMatrix_)));
  #endif
}

void
Window::draw()
{
  glClearColor (1.0, 1.0, 1.0, 0.0); // white
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  setMatrices();

  // this is where the drawing request occurs!
  for (index_t k=0;k<plot_.size();k++)
  {
    plot_[k]->draw();
  }
}

void
Window::write()
{
  // this is where the writing to opengl buffers happens!
  for (index_t k=0;k<plot_.size();k++)
    plot_[k]->write();
}

void
Window::attach( Plot_ptr plot )
{
  // temporary until more generic plots are written
  plot->setWindow(this);
  plot_.push_back( plot );
}

void
Window::run()
{
  // ensure we can capture the escape key being pressed below
  glfwSetInputMode(window_, GLFW_STICKY_KEYS, GL_TRUE);

  // hide the mouse and enable unlimited mouvement
  //glfwSetInputMode(window_, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

  glfwPollEvents();
  glfwSetCursorPos(window_, width_/2, height_/2);

  angles_[0] = 0;
  angles_[1] = 0;

  glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

  // Enable depth test
  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);

  // Cull triangles which normal is not towards the camera
  //glEnable(GL_CULL_FACE);
  glDisable(GL_CULL_FACE);

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  write();

  while (!glfwWindowShouldClose(window_) && glfwGetKey(window_, GLFW_KEY_ESCAPE ) != GLFW_PRESS)
  {
    draw();
    glfwSwapBuffers(window_);
    glfwPollEvents();
  }
}

} // graphics

} // ursa
