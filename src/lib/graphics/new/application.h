#ifndef AVRO_LIB_GRAPHICS_APPLICATION_H_
#define AVRO_LIB_GRAPHICS_APPLICATION_H_

#include "avro_params.h"

#include "graphics/new/gl.h"
#include "graphics/new/vertex_array.h"

#include "graphics/math.h"

#include <memory>
#include <vector>

namespace avro
{

class TopologyBase;

namespace graphics
{

class UniformSet : public ParameterSet {

};

class Shader {
public:
	void setUniform( const std::string& , int u );
	void setUniform( const std::string& , vec3 u );
	void setUniform( const std::string& , mat4 u );
};

class Plot {

public:
	Plot( const TopologyBase& topology ) {

		if (topology.number() == 4) {
			// we need to split up the topology into the boundaries

		}

		// support for parallel meshes too?

	}

	Plot();

	void add( const TopologyBase& topology ) {
		vertex_array_objects_.push_back( std::make_shared<VertexArray>(topology) );
	}

private:
	std::vector< std::shared_ptr<VertexArray> > vertex_array_objects_;
};

class Trackball {
public:

	Trackball( mat4& target ) :
		target_(target)
	{}

	void rotate();
	void translate();

private:
	vec3 center_;
	vec3 radius_;
	mat4& target_;
};

class Camera {

public:
	Camera( const vec3& eye , const vec3& up , const vec3& lookat ,
					real_t fov , index_t width , index_t height );

	vec3& eye() { return eye_; }

private:
	vec3 eye_;
	vec3 up_;
	vec3 lookat_;

	mat4 projection_matrix_;
	mat4 view_matrix_;
};

class Scene {
	/*
	a collection of primitives (points, edges, triangle) that can be rendered
	in a window, or in an HTML canvas
	*/
public:
	typedef std::shared_ptr<Primitive> Primitive_ptr;

	Scene( Manager& manager );

	// create a new primitive
	void add( const TopologyBase& topology );

private:
	Camera camera_;
	std::vector<Primitive_ptr> primitive_;
	Manager& manager_;
};


class Window {
public:

	Window(index_t width , index_t height);
	~Window();

	void render();

	Scene& scene() { return scene_; }
	const Scene& scene() const { return scene_; }

	void mouse_down();
	void mouse_up();
	void mouse_move();
	void key_down();
	void key_up();

private:
	index_t width_;
	index_t height_;
	Scene scene_;
	Interface interface_;
	Trackball trackball_;

	GLFW_Window* window_;
	GLFW_Manager manager_;

	bool picking_;
};

class GLFW_Application : public Application {

public:
	GLFW_Application() :
		Application(manager_)
	{}

	void run();

private:
	GL_Manager manager_;
	Window window_;
};

class Web_Application : public Application {
public:
	Web_Application() :
		Application(manager_)
	{}

	void run();

private:
	WV_Manager manager_;
	Scene scene_;
};

} // graphics

} // avro

#endif
