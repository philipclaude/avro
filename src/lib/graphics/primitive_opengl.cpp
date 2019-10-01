#include "common/tools.h"

#include "graphics/gl.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"

#include "mesh/topology.h"
#include "mesh/vertices.h"

namespace ursa
{

namespace graphics
{

typedef struct
{
  std::vector<float> coordinates;
  std::vector<unsigned int> indices;
} glData;

void
glConvert( const TopologyHolder& topology , glData& data )
{

  index_t nb_elem = topology.nb();
  index_t nb_vert = topology.vertices().nb();
  coord_t dim = topology.vertices().dim();

  const Vertices& vertices = topology.vertices();

  ursa_assert( dim==3 );

  data.coordinates.resize( 3*nb_vert );
  data.indices.resize( 3*nb_elem ); // TODO generalize by using getTriangles...

  index_t n = 0;
  for (index_t k=0;k<nb_elem;k++)
  for (index_t j=0;j<3;j++)
    data.indices[n++] = topology(k,j);

  n = 0;
  for (index_t k=0;k<nb_vert;k++)
  for (index_t j=0;j<dim;j++)
    data.coordinates[n++] = vertices(k,j);
}

void
OpenGLPrimitive::write()
{
  // bind the buffers to the opengl context
  index_t nb_triangles = topology_.nb();
  index_t nb_vertices = topology_.vertices().nb();

  // Create the vertex array object (list of buffers)
  glGenVertexArrays( 1, &vao_ );
  glBindVertexArray(vao_);

  // allocate the buffers
  vbo_.resize(2);
  glGenBuffers(2, vbo_.data() );
  GLuint position = vbo_[0];
  GLuint indices  = vbo_[1];

  // convert the data for the gl
  glData data;
  glConvert( topology_ , data );

  // bind the position buffer
  glBindBuffer(GL_ARRAY_BUFFER, position);
  glBufferData(GL_ARRAY_BUFFER, 3 * nb_vertices * sizeof(float), data.coordinates.data() , GL_STATIC_DRAW);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, (GLubyte *)NULL );

  // bind the index buffer
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indices);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3 * nb_triangles * sizeof(unsigned int), data.indices.data() , GL_STATIC_DRAW);

  // reset the vao bound to the gl
  glBindVertexArray(0);

}

void
OpenGLPrimitive::draw()
{
  ursa_assert( shader_!=NULL );

  shader_->use();

  //shader_->printActiveAttribs();
  //shader_->printActiveUniforms();

  index_t nb_triangles = 2;
  glBindVertexArray(vao_);
  glDrawElements(GL_TRIANGLES, 3 * nb_triangles, GL_UNSIGNED_INT, ((GLubyte *)NULL + (0)));

}

} // graphics

} // ursa
