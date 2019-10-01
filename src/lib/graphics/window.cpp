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

  for (int i=0;i<3;i++)
    mvp_(i,i) = 1;
}

Window::~Window()
{
  glfwDestroyWindow(window_);
}

void
Window::draw()
{
  glClearColor (1.0, 0.0, 0.0, 0.0);

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
  glfwSwapInterval(1);

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
