#include "common/error.h"

#include "graphics/gl.h"
#include "graphics/plot.h"
#include "graphics/plotter.h"
#include "graphics/window.h"

#include "numerics/field.h"

namespace ursa
{

namespace graphics
{

Window::Window( const std::string& title , Plotter* plotter) :
  title_(title),
  plotter_(plotter)
{
  window_ = glfwCreateWindow( width_ , height_ , title_.c_str() , NULL, NULL);
  if (!window_)
  {
    glfwTerminate();
    ursa_assert_not_reached;
  }

  angles_[0] = 0;
  angles_[1] = 0;

  setMatrices();

}

Window::~Window()
{
  glfwDestroyWindow(window_);
}

void
Window::setMatrices()
{
  // glfwGetTime is called only once, the first time this function is called
	static double lastTime = glfwGetTime();


  	// Compute time difference between current and last frame
    #if 0
  	double currentTime = glfwGetTime();
  	float deltaTime = float(currentTime - lastTime);
    deltaTime = 0.f;

  	// Get mouse position
  	double xpos, ypos;
  	glfwGetCursorPos(window_, &xpos, &ypos);

  	// Reset mouse position for next frame
  	//glfwSetCursorPos(window_, 1024/2, 768/2);

  	// Compute new orientation
  	angles_[0] += mouseSpeed_ * float(width_/2  - xpos );
  	angles_[1] += mouseSpeed_ * float(height_/2 - ypos );

  	// Direction : Spherical coordinates to Cartesian coordinates conversion
  	glm::vec3 direction(
  		cos(angles_[1]) * sin(angles_[0]),
  		sin(angles_[1]),
  		cos(angles_[1]) * cos(angles_[0])
  	);

  	// Right vector
  	glm::vec3 right = glm::vec3(
  		sin(angles_[0] - 3.14f/2.0f),
  		0,
  		cos(angles_[0] - 3.14f/2.0f)
  	);

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

  	// Projection matrix : 45° Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
  	projMatrix_ = glm::perspective(glm::radians(fov_), 4.0f / 3.0f, 0.1f, 100.0f);

  	// Camera matrix
  	viewMatrix_ = glm::lookAt(
  								position_,           // Camera is here
  								position_+direction, // and looks here : at the same position, plus "direction"
  								up                  // Head is up (set to 0,-1,0 to look upside-down)
  						   );

   #else
   projMatrix_ = glm::perspective( glm::radians(fov_) , (float) width_  / height_ , 0.1f , 40.0f );
   projMatrix_ = glm::translate( projMatrix_ , vec3(0,0,-1) );

   #endif

    modelMatrix_ = mat4(1.0);

    mvp_ = projMatrix_*viewMatrix_*modelMatrix_;

  	// For the next frame, the "last time" will be "now"
  	//lastTime = currentTime;
}

void
Window::draw()
{
  glClearColor (1.0, 1.0, 1.0, 0.0); // white

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
  glfwMakeContextCurrent(window_);
  //glfwSwapInterval(1);

  glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

  // Enable depth test
  glEnable(GL_DEPTH_TEST);
  // Accept fragment if it closer to the camera than the former one
  glDepthFunc(GL_LESS);

  // Cull triangles which normal is not towards the camera
  glEnable(GL_CULL_FACE);

  write();

  while (!glfwWindowShouldClose(window_))
  {

    draw(); // callback to window virtual function
    /*mat4x4 m, p, mvp;
    mat4x4_identity(m);
    mat4x4_rotate_Z(m, m, (float) glfwGetTime());
    mat4x4_ortho(p, -ratio, ratio, -1.f, 1.f, 1.f, -1.f);
    mat4x4_mul(mvp, p, m);

    glUseProgram(program);
    glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (const GLfloat*) &mvp);
    glBindVertexArray(vertex_array);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    */
    glfwSwapBuffers(window_);
    glfwPollEvents();
  }
}

} // graphics

} // ursa
