#include "graphics/application.h"

namespace avro
{

namespace graphics
{

void
ApplicationBase::write()
{
  // write all the data present in the scene graph
  for (index_t k=0;k<scene_.size();k++)
    scene_[k].write(manager_);
}

/*
void
OpenGL_Manager::draw()
{
  // draws with the current shader (selected by the scene graph)
  std::map<index_t,index_t>::iterator it;

  for (it=vao_triangles_.begin();it!=vao_triangles_.end();it++)
  {
    index_t vao = it->first;
    index_t nb = it->second;
    bool active = vao_triangles_active_[vao];
    ShaderProgram* shader = vao_triangles_shader_[vao];

    if (!active) continue;

    shader->use();
    GL_CALL( glBindVertexArray(vao) );
    GL_CALL( glDrawElements(GL_TRIANGLES,nb, GL_UNSIGNED_INT , 0 ) )
  }

  for (it=vao_edges_.begin();it!=vao_edges_.end();it++)
  {
    index_t vao = it->first;
    index_t nb = it->second;
    bool active = vao_edges_active_[vao];
    ShaderProgram* shader = vao_edges_shader_[vao];

    if (!active) continue;

    shader->use();
    GL_CALL( glBindVertexArray(vao) );
    GL_CALL( glDrawElements( GL_LINES , nb , GL_UNSIGNED_INT , 0 ) )
  }

  for (it=vao_points_.begin();it!=vao_points_.end();it++)
  {
    index_t vao = it->first;
    index_t nb = it->second;
    bool active = vao_points_active_[vao];
    ShaderProgram* shader = vao_points_shader_[vao];

    if (!active) continue;

    shader->use();
    GL_CALL( glBindVertexArray(vao) );
    GL_CALL( glPointSize(10.0f) );
    GL_CALL( glDrawArrays( GL_POINTS , 0 , nb ) );
  }
}
*/

void
ShaderLibrary::set_matrices( SceneGraph& scene )
{
  // go through all the active shaders and assign the MVP and normalMatrix
  std::map<std::string,ShaderProgram>::iterator it;
  for (it=shaders_.begin();it!=shaders_.end();it++)
  {
    ShaderProgram& shader = it->second;
    shader.setUniform("MVP" , scene.mvp_matrix() );
    shader.setUniform("u_normalMatrix" , scene.normal_matrix() );
  }
}

void
OpenGL_Manager::draw( Primitive& primitive , TransformFeedbackResult* feedback )
{
  // indicate to the gl that we want to use the shader
  ShaderProgram shader = *shader_.at(&primitive);
  shader.use();

  if (feedback!=nullptr)
  {
    GLuint& query = feedback->query();
    GLuint& buffer = feedback->buffer();

    glGenQueries(1, &query );
    GL_CALL( glBeginQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN , query ) );

    GL_CALL( glGenBuffers( 1 , &buffer ) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , buffer ) );
    GL_CALL( glBufferData( GL_TRANSFORM_FEEDBACK_BUFFER , 8*primitive.nb_triangles()*sizeof(GLfloat) , NULL , GL_DYNAMIC_COPY ) );
    GL_CALL( glBindBufferBase( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , buffer ) );
    GL_CALL( glBindVertexArray( feedback->vao() ) );

    GL_CALL( glEnable(GL_RASTERIZER_DISCARD) );
    GL_CALL( glBeginTransformFeedback(GL_TRIANGLES) );
  }

  // bind the vao associated with this primitive
  if (primitive.number()>=2)
  {
    // draw the triangles
    if (primitive.triangles_on())
    {
      GL_CALL( glBindVertexArray( vao_triangles_.at(&primitive) ) );
      GL_CALL( glDrawElements(GL_TRIANGLES, 3*primitive.nb_triangles() , GL_UNSIGNED_INT , 0 ) )
    }
  }

  if (feedback!=nullptr)
  {
    GL_CALL( glEndTransformFeedback() );
    GL_CALL( glDisable(GL_RASTERIZER_DISCARD) );

    GLuint nb_primitives = 0;
    GL_CALL( glEndQuery( GL_TRANSFORM_FEEDBACK_PRIMITIVES_WRITTEN ) );
    GL_CALL( glGetQueryObjectuiv( feedback->query() , GL_QUERY_RESULT , &nb_primitives ) );
    printf("primitives written = %u\n",nb_primitives);

    // there are 4 coordinates for every output and 3 vertices per prim (color + coord)
    index_t nb_data_per_prim = 4*3*2;
    index_t size = nb_primitives*nb_data_per_prim;

    GLfloat* feedbackBuffer = (GLfloat*) malloc( size*sizeof(GLfloat) );
    GL_CALL( glBindBuffer( GL_TRANSFORM_FEEDBACK_BUFFER , feedback->buffer() ) );
    GL_CALL( glGetBufferSubData( GL_TRANSFORM_FEEDBACK_BUFFER , 0 , size*sizeof(GLfloat) , feedbackBuffer ) );

    //printBuffer(feedbackBuffer, nb_primitives , nb_data_per_prim );
  }

  if (primitive.number()>=1)
  {
    if (primitive.edges_on())
    {
      GL_CALL( glBindVertexArray( vao_edges_.at(&primitive) ) );
      GL_CALL( glDrawElements( GL_LINES , primitive.nb_edges()*2 , GL_UNSIGNED_INT , 0 ) )
    }
  }

  // draw the points
  if (primitive.points_on())
  {
    GL_CALL( glBindVertexArray(vao_points_.at(&primitive) ) );
    GL_CALL( glPointSize(10.0f) );
    GL_CALL( glDrawArrays( GL_POINTS , 0 , primitive.nb_points() ) );
  }

  // reset the vao bound to the gl
  GL_CALL( glBindVertexArray(0) );
}

void
OpenGL_Manager::draw( SceneGraph& scene , TransformFeedbackResult* feedback )
{
  if (!scene.update()) return;

  shaders_.set_matrices(scene);

  for (index_t k=0;k<scene.nb_primitives();k++)
  {
    draw( scene.primitive(k) , feedback );
  }

  scene.set_update(false);
}

template class Application<GLFW_Interface<OpenGL_Manager>>;
template class Application<GLFW_Interface<Vulkan_Manager>>;

} // graphics

} // avro
