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
  data.normals.resize( 3*nb_vert , 0 );

  index_t n = 0;
  for (index_t k=0;k<nb_elem;k++)
  for (index_t j=0;j<nv;j++)
    data.indices[n++] = topology_(k,j);

  n = 0;
  for (index_t k=0;k<nb_vert;k++)
  {
    data.colours[3*k  ] = 255;
    data.colours[3*k+1] = 255;
    data.colours[3*k+2] = 0;

    for (index_t j=0;j<dim;j++)
    {
      data.coordinates[n] = vertices(k,j);
      n++;
    }
  }

  // compute the normals
  if (topology_.number()==2)
  {
    real_t x1,y1,z1,x2,y2,z2,x3,y3,z3;
    real_t dis;
    real_t xnor,ynor,znor;
    std::vector<index_t> count( nb_vert , 0 );
    for (index_t k=0;k<nb_elem;k++)
    {
      x1 = vertices( topology_(k,0) , 0 );
      y1 = vertices( topology_(k,0) , 1 );
      z1 = vertices( topology_(k,0) , 2 );
      x2 = vertices( topology_(k,1) , 0 );
      y2 = vertices( topology_(k,1) , 1 );
      z2 = vertices( topology_(k,1) , 2 );
      x3 = vertices( topology_(k,2) , 0 );
      y3 = vertices( topology_(k,2) , 1 );
      z3 = vertices( topology_(k,2) , 2 );

      xnor = (y3-y1)*(z2-z1)-(y2-y1)*(z3-z1);
      ynor = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1);
      znor = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);

      dis  = sqrtf(xnor*xnor + ynor*ynor + znor*znor);

      if (dis != 0.0f)
      {
        data.normals[ 3*topology_(k,0) + 0 ] += xnor/dis;
        data.normals[ 3*topology_(k,0) + 1 ] += ynor/dis;
        data.normals[ 3*topology_(k,0) + 2 ] += znor/dis;
        data.normals[ 3*topology_(k,1) + 0 ] += xnor/dis;
        data.normals[ 3*topology_(k,1) + 1 ] += ynor/dis;
        data.normals[ 3*topology_(k,1) + 2 ] += znor/dis;
        data.normals[ 3*topology_(k,2) + 0 ] += xnor/dis;
        data.normals[ 3*topology_(k,2) + 1 ] += ynor/dis;
        data.normals[ 3*topology_(k,2) + 2 ] += znor/dis;

        count[ topology_(k,0) ]++;
        count[ topology_(k,1) ]++;
        count[ topology_(k,2) ]++;
      }
    }

    // normalize again
    for (index_t k=0;k<nb_vert;k++)
    {
      if (count[k] <= 1) continue;
      dis  = count[k];
      xnor = data.normals[3*k  ] / dis;
      ynor = data.normals[3*k+1] / dis;
      znor = data.normals[3*k+2] / dis;
      dis  = sqrtf(xnor*xnor + ynor*ynor + znor*znor);
      data.normals[3*k  ] = xnor/dis;
      data.normals[3*k+1] = ynor/dis;
      data.normals[3*k+2] = znor/dis;
    }
  }
  else if (topology_.number()==1)
  {
    printf("setting color for lines!!\n");
    data.colours.resize( 3*nb_vert , 255.0f );
  }

  // TODO add extra stuff stored in the topology
  const Fields& fields = topology_.fields();
  if (fields.has("vertex_normals"))
  {
    printf("add vertex normals!!!\n");
  }
  if (fields.has("vertex_uv"))
  {
    printf("add vertex uv values!!\n");
  }

  if (fields.has("triangle_normals"))
  {
    printf("add triangle normals!!!\n");
  }
  if (fields.has("triangle_uv"))
  {
    printf("add triangle uv values!!\n");
  }
}

void
OpenGLPrimitive::write()
{
  // bind the buffers to the opengl context
  index_t nb_elem = topology_.nb();
  index_t nb_vertices = topology_.vertices().nb();
  coord_t nv = topology_.number()+1;

  // allocate the buffers
  vbo_.resize(4);
  GL_CALL( glGenBuffers( vbo_.size() , vbo_.data() ) );
  GLuint& indices  = vbo_[0];
  GLuint& position = vbo_[1];
  GLuint& colour   = vbo_[2];
  GLuint& normal   = vbo_[3];

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

  // bind the normal buffer
  if (topology_.number()==2)
  {
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, normal) );
    GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_vertices * sizeof(GLfloat), data_.normals.data() , GL_STATIC_DRAW) );
  }

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

  if (topology_.number()==2)
  {
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, normal ) );
    GL_CALL( glVertexAttribPointer( 2, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(2) );
  }

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );

}

void
OpenGLPrimitive::draw()
{
  ursa_assert( shader_!=NULL );

  index_t nb_elem = topology_.nb();
  index_t nv = topology_.number()+1;

  // indicate to the gl that we want to use the shader
  shader_->use();

  // assign the uniforms for the shader programs
  shader_->setUniform("MVP" , window_->mvp() );
  shader_->setUniform("u_normalMatrix" , window_->normal() );

  // bind the vao associated with this primitive
  GL_CALL( glBindVertexArray(vao_) );
  if (topology_.number()==1 && visible_)
  {
    shader_->setUniform("wLight" , 0.0f );
    shader_->setUniform("wColor" , 0.0f );
    GL_CALL( glDrawElements( GL_LINES , nb_elem*nv , GL_UNSIGNED_INT , 0 ) )
  }
  else if (topology_.number()==2 && visible_)
  {
    GL_CALL( glDrawElements(GL_TRIANGLES,nb_elem*nv, GL_UNSIGNED_INT , 0 ) )
  }

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

} // graphics

} // ursa
