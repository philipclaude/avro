#ifndef AVRO_LIB_GRAPHICS_MANAGERS_H_
#define AVRO_LIB_GRAPHICS_MANAGERS_H_

#include "graphics/gl.h"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace avro
{

namespace graphics
{

class VertexAttributeObject;
class ShaderProgram;
class Shaders;

class OpenGL4_Manager {
public:
	OpenGL4_Manager();
	~OpenGL4_Manager();

	void write( VertexAttributeObject& vao );

	void track_buffer(gl_index b) { buffers_.push_back(b); }
	void track_texture(gl_index t) { textures_.push_back(t); }

private:
	// keep track of things to delete
	std::vector<gl_index> buffers_;
	std::vector<gl_index> textures_;
	std::vector<gl_index> vao_;

  //std::map< std::string , std::shared_ptr<ShaderProgram> > shaders_;

public:
	static std::shared_ptr<Shaders> shaders;
};


} // graphics

} // avro

#endif
