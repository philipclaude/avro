#ifndef AVRO_LIB_GRAPHICS_MANAGERS_H_
#define AVRO_LIB_GRAPHICS_MANAGERS_H_

#include "graphics/gl.h"

#include <map>
#include <string>

namespace avro
{

namespace graphics
{

class VertexAttributeObject;
class ShaderProgram;

class OpenGL4_Manager {
public:
	OpenGL4_Manager();
	~OpenGL4_Manager();

	void write( VertexAttributeObject& vao );

private:
	// keep track of things to delete
	std::vector<gl_index> buffers_;
	std::vector<gl_index> textures_;
	std::vector<gl_index> vao_;

  std::map< std::string , std::shared_ptr<ShaderProgram> > shaders_;
};


} // graphics

} // avro

#endif