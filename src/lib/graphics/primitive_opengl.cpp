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
OpenGLPrimitive::convert()
{
  index_t n = 0;
  index_t nb_vertices = topology_.vertices().nb();
  index_t nb_triangles = 0;
  index_t nb_edges = 0;
  coord_t dim = topology_.vertices().dim();
  coord_t nv = topology_.number()+1;
  coord_t number = topology_.number();

  const Vertices& vertices = topology_.vertices();

  ursa_assert( dim==3 );

  // allocate the vertex data
  points_.resize( 3*nb_vertices);
  colors_.resize( 3*nb_vertices, 0 );
  normals_.resize( 3*nb_vertices, 0 );

  if (number>=1)
  {
    // get the edges from the topology
    std::vector<index_t> edges;
    topology_.getEdges( edges );
    nb_edges = edges_.size()/2;
    edges_.assign( edges.begin() , edges.end() );
  }
  if (number>=2)
  {
    // get the triangles from the topology (TODO)
    ursa_assert( number==2 );
    nb_triangles = topology_.nb();
    triangles_.resize( 3*nb_triangles );

    // save the triangles
    n = 0;
    for (index_t k=0;k<nb_triangles;k++)
    for (index_t j=0;j<3;j++)
      triangles_[n++] = topology_(k,j);
  }

  // save the vertex data that gets passed to the shaders
  n = 0;
  for (index_t k=0;k<nb_vertices;k++)
  {
    colors_[3*k  ] = 255;
    colors_[3*k+1] = 255;
    colors_[3*k+2] = 0;

    for (index_t j=0;j<dim;j++)
    {
      points_[n] = vertices(k,j);
      n++;
    }
  }

  // compute the normals
  if (topology_.number()==2)
  {
    real_t x1,y1,z1,x2,y2,z2,x3,y3,z3;
    real_t dis;
    real_t xnor,ynor,znor;
    std::vector<index_t> count( nb_vertices, 0 );
    for (index_t k=0;k<nb_triangles;k++)
    {
      x1 = vertices( triangles_[3*k  ] , 0 );
      y1 = vertices( triangles_[3*k  ] , 1 );
      z1 = vertices( triangles_[3*k  ] , 2 );
      x2 = vertices( triangles_[3*k+1] , 0 );
      y2 = vertices( triangles_[3*k+1] , 1 );
      z2 = vertices( triangles_[3*k+1] , 2 );
      x3 = vertices( triangles_[3*k+2] , 0 );
      y3 = vertices( triangles_[3*k+2] , 1 );
      z3 = vertices( triangles_[3*k+2] , 2 );

      xnor = (y3-y1)*(z2-z1)-(y2-y1)*(z3-z1);
      ynor = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1);
      znor = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);

      dis  = sqrtf(xnor*xnor + ynor*ynor + znor*znor);

      if (dis != 0.0f)
      {
        normals_[ 3*triangles_[3*k  ] + 0 ] += xnor/dis;
        normals_[ 3*triangles_[3*k  ] + 1 ] += ynor/dis;
        normals_[ 3*triangles_[3*k  ] + 2 ] += znor/dis;
        normals_[ 3*triangles_[3*k+1] + 0 ] += xnor/dis;
        normals_[ 3*triangles_[3*k+1] + 1 ] += ynor/dis;
        normals_[ 3*triangles_[3*k+1] + 2 ] += znor/dis;
        normals_[ 3*triangles_[3*k+2] + 0 ] += xnor/dis;
        normals_[ 3*triangles_[3*k+2] + 1 ] += ynor/dis;
        normals_[ 3*triangles_[3*k+2] + 2 ] += znor/dis;

        count[ triangles_[3*k  ] ]++;
        count[ triangles_[3*k+1] ]++;
        count[ triangles_[3*k+2] ]++;
      }
    }

    // normalize again
    for (index_t k=0;k<nb_vertices;k++)
    {
      if (count[k] <= 1) continue;
      dis  = count[k];
      xnor = normals_[3*k  ] / dis;
      ynor = normals_[3*k+1] / dis;
      znor = normals_[3*k+2] / dis;
      dis  = sqrtf(xnor*xnor + ynor*ynor + znor*znor);
      normals_[3*k  ] = xnor/dis;
      normals_[3*k+1] = ynor/dis;
      normals_[3*k+2] = znor/dis;
    }
  }
  else if (topology_.number()==1)
  {
    printf("setting color for lines!!\n");
    colors_.resize( 3*nb_vertices, 255.0f );
  }

  // TODO add extra stuff stored in the topology
  /*
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
  */
}

void
OpenGLPrimitive::write()
{
  // bind the buffers to the opengl context
  index_t nb_vertices = topology_.vertices().nb();

  // allocate the buffers
  std::vector<GLuint> vbo(6);
  GL_CALL( glGenBuffers( vbo.size() , vbo.data() ) );
  GLuint& position  = vbo[0];
  GLuint& colour    = vbo[1];
  GLuint& normal    = vbo[2];
  GLuint& triangles = vbo[3];
  GLuint& edges     = vbo[4];
  GLuint& points    = vbo[5];

  // convert the data for the gl
  convert();

  if (topology_.number()>=2)
  {
    // bind the triangles
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangles) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles_.size() * sizeof(GLuint), triangles_.data() , GL_STATIC_DRAW) );
  }

  // bind the position buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, position) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_vertices * sizeof(GLfloat), points_.data() , GL_STATIC_DRAW) );

  // bind the colour buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, colour) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_vertices * sizeof(GLfloat), colors_.data() , GL_STATIC_DRAW) );

  // bind the normal buffer
  if (topology_.number()>=2)
  {
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, normal) );
    GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_vertices * sizeof(GLfloat), normals_.data() , GL_STATIC_DRAW) );
  }

  // bind the triangle data
  GL_CALL( glGenVertexArrays( 1, &vao_triangles_ ) );
  GL_CALL( glBindVertexArray(vao_triangles_) );

  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, triangles ) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, colour ) );
  GL_CALL( glVertexAttribPointer( 1, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(1) );

  if (topology_.number()>=2)
  {
    GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, normal ) );
    GL_CALL( glVertexAttribPointer( 2, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
    GL_CALL( glEnableVertexAttribArray(2) );
  }

  // bind the edges
  GL_CALL( glGenVertexArrays( 1, &vao_edges_ ) );
  GL_CALL( glBindVertexArray(vao_edges_) );

  GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, edges ) );
  GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, edges_.size() * sizeof(GLuint), edges_.data() , GL_STATIC_DRAW) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

  // bind the points
  GL_CALL( glGenVertexArrays( 1, &vao_points_ ) );
  GL_CALL( glBindVertexArray(vao_points_) );

  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, position ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );

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
  if (topology_.number()>=2)
  {
    // draw the triangles
    if (triangles_on_)
    {
      GL_CALL( glBindVertexArray(vao_triangles_) );
      GL_CALL( glDrawElements(GL_TRIANGLES,triangles_.size(), GL_UNSIGNED_INT , 0 ) )
    }
  }
  if (topology_.number()>=1)
  {
    if (edges_on_)
    {
      GL_CALL( glBindVertexArray(vao_edges_) );
      GL_CALL( glDrawElements( GL_LINES , edges_.size() , GL_UNSIGNED_INT , 0 ) )
    }
  }
  if (topology_.number()>=0)
  {
    // draw the points
    if (points_on_)
    {
      GL_CALL( glBindVertexArray(vao_points_) );
      GL_CALL( glPointSize(10.0f) );
      GL_CALL( glDrawArrays( GL_POINTS , 0 , points_.size()/3 ) );
    }

  }

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

} // graphics

} // ursa
