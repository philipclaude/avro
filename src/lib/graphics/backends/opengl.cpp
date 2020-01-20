#include "common/tools.h"

#include "graphics/colormap.h"
#include "graphics/gl.h"
#include "graphics/primitive.h"
#include "graphics/shader.h"
#include "graphics/window.h"

#include "library/eps.h"

#include "mesh/decomposition.h"
#include "mesh/field.hpp"
#include "mesh/points.h"
#include "mesh/topology.h"

namespace avro
{

namespace graphics
{

void
Primitive::convert()
{
  index_t nb_triangles = 0;
  coord_t dim0 = topology_.points().dim();
  coord_t number = topology_.number();

  // get the triangles from the topology
  std::shared_ptr<SimplicialDecompositionBase> pdecomposition;
  if (topology_.type_name()=="simplex")
    pdecomposition = std::make_shared<SimplicialDecomposition<Simplex>>(*dynamic_cast<const Topology<Simplex>*>(&topology_));
  else
    pdecomposition = std::make_shared<SimplicialDecomposition<Polytope>>(*dynamic_cast<const Topology<Polytope>*>(&topology_));
  pdecomposition->extract();
  const SimplicialDecompositionBase& decomposition = *pdecomposition.get();

  const Points& points = decomposition.points();
  std::vector<index_t> triangles;
  std::vector<index_t> parents;
  decomposition.get_simplices(2,triangles,parents); // setting 2 retrieves triangles
  index_t nb_points = triangles.size(); // duplicated points for cell-based colors

  // allocate the vertex data
  colors_.resize( 3*nb_points, 0 );
  normals_.resize( 3*nb_points, 0 );

  std::vector<index_t> point_map0( points.nb() );
  std::vector<index_t> point_map1;
  for (index_t k=0;k<triangles.size();k++)
  {
    for (index_t j=0;j<dim0;j++)
      points_.push_back( points[triangles[k]][j] );
    for (index_t j=dim0;j<3;j++)
      points_.push_back( 0.0 );

    point_map0[ triangles[k] ] = triangles_.size();
    point_map1.push_back( triangles[k] );

    if (number>=2)
      triangles_.push_back( triangles_.size() );
  }
  nb_triangles = triangles_.size()/3;

  if (number>=1)
  {
    // get the edges from the topology
    std::vector<index_t> edges;
    topology_.get_edges( edges );
    for (index_t k=0;k<edges.size();k++)
      edges_.push_back( point_map0[ edges[k] ] );
  }

  std::vector<real_t> U;
  real_t umin,umax;

  bool constant_color = false;
  float color[3];
  const Fields& fields = topology_.fields();
  if (fields.has(active_))
  {
    // retrieve the active field
    umin = fields[active_].min(rank_);
    umax = fields[active_].max(rank_);

    // evaluate the active field (with rank) on the triangulation points
    const std::vector<index_t>& point_parents = decomposition.point_parents();
    const Table<real_t>& reference_coordinates = decomposition.reference_coordinates();
    fields[active_].evaluate( rank_ , point_parents , reference_coordinates , U );
  }
  else
  {
    // constant color
    umin = 0.0;
    umax = 1.0;
    U.resize( decomposition.points().nb() , 0.5 );
    color[0] = 0.7;
    color[1] = 0.7;
    color[2] = 0.7;
    constant_color = true;
  }
  float lims[2] = {float(umin),float(umax)};

  // compute the color of the (unduplicated) triangulation points
  Colormap colormap;
  colormap.set_limits(lims);
  std::vector<index_t> color0( points.nb()*3 );
  for (index_t k=0;k<triangles.size()/3;k++)
  {
    // retrieve the parent
    for (index_t j=0;j<3;j++)
    {
      index_t p = triangles[3*k+j];
      real_t  u = U[p];
      if (!constant_color)
        colormap.map( u , color );

      color0[3*p  ] = 255*color[0];
      color0[3*p+1] = 255*color[1];
      color0[3*p+2] = 255*color[2];
    }
  }

  // map the colors to all the vertices
  for (index_t k=0;k<nb_points;k++)
  {
    for (index_t j=0;j<3;j++)
      colors_[3*k+j] = color0[3*point_map1[k]+j];
  }

  // compute the normals
  if (topology_.number()>=2)
  {
    real_t x1,y1,z1,x2,y2,z2,x3,y3,z3;
    real_t dis;
    real_t xnor,ynor,znor;
    std::vector<index_t> count( nb_points, 0 );
    for (index_t k=0;k<nb_triangles;k++)
    {
      x1 = points_[ 3*triangles_[3*k  ] + 0 ];
      y1 = points_[ 3*triangles_[3*k  ] + 1 ];
      z1 = points_[ 3*triangles_[3*k  ] + 2 ];
      x2 = points_[ 3*triangles_[3*k+1] + 0 ];
      y2 = points_[ 3*triangles_[3*k+1] + 1 ];
      z2 = points_[ 3*triangles_[3*k+1] + 2 ];
      x3 = points_[ 3*triangles_[3*k+2] + 0 ];
      y3 = points_[ 3*triangles_[3*k+2] + 1 ];
      z3 = points_[ 3*triangles_[3*k+2] + 2 ];

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
    for (index_t k=0;k<nb_points;k++)
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
    colors_.resize( 3*nb_points, 255.0f );
  }
}

void
OpenGLPrimitive::write()
{
  printf("writing!!\n");
  
  // bind the buffers to the opengl context
  transform_feedback_ = false;

  // allocate the buffers
  std::vector<GLuint> vbo(7);
  GL_CALL( glGenBuffers( vbo.size() , vbo.data() ) );
  GLuint& position  = vbo[0];
  GLuint& colour    = vbo[1];
  GLuint& normal    = vbo[2];
  GLuint& triangles = vbo[3];
  GLuint& edges     = vbo[4];
  GLuint& points    = vbo[5]; UNUSED(points);
  GLuint& feedback  = vbo[6];

  // convert the data for the gl
  convert();

  index_t nb_points = points_.size()/3;

  if (topology_.number()>=2)
  {
    // bind the triangles
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, triangles) );
    GL_CALL( glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles_.size() * sizeof(GLuint), triangles_.data() , GL_STATIC_DRAW) );
  }

  // bind the position buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, position) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_points * sizeof(GLfloat), points_.data() , GL_STATIC_DRAW) );

  // NEW
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, feedback) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_points * sizeof(GLfloat) , NULL , GL_STATIC_COPY) );
  glBindBuffer(GL_ARRAY_BUFFER,0);
  // END NEW

  // bind the colour buffer
  GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, colour) );
  GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_points * sizeof(GLfloat), colors_.data() , GL_STATIC_DRAW) );

  // bind the normal buffer
  if (topology_.number()>=2)
  {
    GL_CALL( glBindBuffer(GL_ARRAY_BUFFER, normal) );
    GL_CALL( glBufferData(GL_ARRAY_BUFFER, 3 * nb_points * sizeof(GLfloat), normals_.data() , GL_STATIC_DRAW) );
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

  // bind the feedback buffer
  // NEW
  GL_CALL( glGenVertexArrays( 1, &vao_feedback_ ) );
  GL_CALL( glBindVertexArray( vao_feedback_ ) );
  GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , feedback ) );
  GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) );
  GL_CALL( glEnableVertexAttribArray(0) );
  // END NEW

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

void
printBuffer( GLfloat* buffer , index_t nb_prim , index_t nb_data_per_prim  )
{
  index_t n = 0;
  std::vector<real_t> triangles;
  std::vector<real_t> colors;
  for (index_t k=0;k<nb_prim;k++)
  {
    //printf("primitve %lu\n",k);
    for (index_t j=0;j<3;j++)
    {
      for (coord_t d=0;d<3;d++)
        triangles.push_back( buffer[n+d] );
      for (coord_t d=0;d<3;d++)
        colors.push_back( buffer[n+4+d] );

      /*printf("v = (%g,%g,%g,%g), c = (%g,%g,%g,%g)\n",
              buffer[n  ],buffer[n+1],buffer[n+2],buffer[n+3] ,
              buffer[n+4],buffer[n+5],buffer[n+6],buffer[n+7]);*/
      n += 8;
    }
  }

  library::epsFile eps;
  int viewport[4] = {0,0,1024,640};
  eps.set_viewport(viewport);
  eps.add_triangles( triangles , colors );
  eps.write( "test.eps" );
}

void
OpenGLPrimitive::draw()
{
  avro_assert( shader_!=NULL );

  // indicate to the gl that we want to use the shader
  shader_->use();

  // assign the uniforms for the shader programs
  shader_->setUniform("MVP" , window_->mvp() );
  shader_->setUniform("u_normalMatrix" , window_->normal() );

  GLuint query;
  GLuint buffer;
  if (transform_feedback_)
  {
    glGenQueries(1, &query);
    GL_CALL( glBeginQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN , query ) );

    GL_CALL( glGenBuffers( 1 , &buffer ) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glBufferData( GL_TRANSFORM_FEEDBACK_BUFFER , 24*triangles_.size()*sizeof(GLfloat) , NULL , GL_DYNAMIC_COPY ) );
    GL_CALL( glBindBufferBase( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , buffer ) );
    GL_CALL( glBindVertexArray( vao_feedback_ ) );

    GL_CALL( glEnable(GL_RASTERIZER_DISCARD) );
    GL_CALL( glBeginTransformFeedback(GL_TRIANGLES) );
  }

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

  if (transform_feedback_)
  {
    GL_CALL( glEndTransformFeedback() );
    GL_CALL( glDisable(GL_RASTERIZER_DISCARD) );

    GLuint nb_primitives = 0;
    GL_CALL( glEndQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN ) );
    GL_CALL( glGetQueryObjectuiv( query , GL_QUERY_RESULT , &nb_primitives ) );
    printf("primitives written = %u\n",nb_primitives);

    // there are 4 coordinates for every output and 3 vertices per prim (color + coord)
    index_t nb_data_per_prim = 4*3*2;
    index_t size = nb_primitives*nb_data_per_prim;

    GLfloat* feedbackBuffer = (GLfloat*) malloc( size*sizeof(GLfloat) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glGetBufferSubData( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , size*sizeof(GLfloat) , feedbackBuffer ) );

    printBuffer(feedbackBuffer, nb_primitives , nb_data_per_prim );
  }

  if (topology_.number()>=1)
  {
    if (edges_on_)
    {
      GL_CALL( glBindVertexArray(vao_edges_) );
      GL_CALL( glDrawElements( GL_LINES , edges_.size() , GL_UNSIGNED_INT , 0 ) )
    }
  }
  // draw the points
  if (points_on_)
  {
    GL_CALL( glBindVertexArray(vao_points_) );
    GL_CALL( glPointSize(10.0f) );
    GL_CALL( glDrawArrays( GL_POINTS , 0 , points_.size()/3 ) );
  }

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

} // graphics

} // avro
