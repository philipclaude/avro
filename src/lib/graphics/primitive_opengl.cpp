#include "common/tools.h"

#include "graphics/gl.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"
#include "graphics/window.h"

#include "mesh/topology.h"
#include "mesh/vertices.h"

namespace ursa
{

namespace graphics
{

void
OpenGLPrimitive::convert( glData& data )
{

  index_t nb_elem = topology_.nb();
  index_t nb_vert = topology_.vertices().nb();
  coord_t dim = topology_.vertices().dim();
  coord_t nv = topology_.number()+1;

  const Vertices& vertices = topology_.vertices();

  ursa_assert( dim==3 );

  data.coordinates.resize( 3*nb_vert );
  data.indices.resize( nv*nb_elem ); // TODO generalize by using getTriangles...
  data.colours.resize( 3*nb_vert, 0 );

  index_t n = 0;
  for (index_t k=0;k<nb_elem;k++)
  for (index_t j=0;j<nv;j++)
    data.indices[n++] = topology_(k,j);

  n = 0;
  for (index_t k=0;k<nb_vert;k++)
  {
    data.colours[3*k] = k*1.0f/nb_vert;
    for (index_t j=0;j<dim;j++)
    {
      //data.colours[n] = (float) (k+1)*(j+1) / (3*nb_vert);
      data.colours[n] = vertices(k,0);
      data.coordinates[n++] = vertices(k,j);
    }
  }
  printInline(data.colours);
}

void
OpenGLPrimitive::write()
{
  glViewport(0,0,1,1);

  // bind the buffers to the opengl context
  index_t nb_elem = topology_.nb();
  index_t nb_vertices = topology_.vertices().nb();
  coord_t nv = topology_.number()+1;

  printf("nb_elem = %lu\n",nb_elem);

  // allocate the buffers
  vbo_.resize(3);
  GL_CALL( glGenBuffers( vbo_.size() , vbo_.data() ) );
  GLuint& indices  = vbo_[0];
  GLuint& position = vbo_[1];
  GLuint& colour   = vbo_[2];

  // convert the data for the gl
  convert( data_ );

  ursa_assert( data_.indices.size()==nv*nb_elem );

  // bind the index buffer
  GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, nv * nb_elem * sizeof(GLuint), data_.indices.data() , GL_STATIC_DRAW) );

  // bind the position buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, position) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_vertices * sizeof(GLfloat), data_.coordinates.data() , GL_STATIC_DRAW) );

  // bind the colour buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, colour) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_vertices * sizeof(GLfloat), data_.colours.data() , GL_STATIC_DRAW) );

  // Create the vertex array object (list of buffers)
  GL_CALL( glGenVertexArrays( 1, &vao_ ) );
  GL_CALL( glBindVertexArray(vao_) );

  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, indices ) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, colour ) );
  GL_CALL( glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(1) );

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );

}

void
OpenGLPrimitive::draw()
{
  ursa_assert( shader_!=NULL );

  index_t nb_elem = topology_.nb();
  index_t nv = topology_.number()+1;
  index_t nb_vert = topology_.vertices().nb();
  coord_t dim = 3;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  shader_->use();

  shader_->setUniform("MVP" , window_->mvp() );

  //shader_->printActiveAttribs();
  //shader_->printActiveUniforms();

  //shader_->setUniform("MVP",window_->mvp());

  //std::cout << window_->mvp() << std::endl;

  //GL_CALL( glClear(GL_COLOR_BUFFER_BIT) );

  GL_CALL( glBindVertexArray(vao_) );
  if (topology_.number()==1)
    GL_CALL( glDrawElements( GL_LINES , nb_elem*nv , GL_UNSIGNED_INT , 0 ) )
  else if (topology_.number()==2)
    GL_CALL( glDrawElements(GL_TRIANGLES,nb_elem*nv, GL_UNSIGNED_INT , 0 ) )
  else
    ursa_implement;
  GL_CALL( glBindVertexArray(0) );

}

} // graphics

} // ursa
