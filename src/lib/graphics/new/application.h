#ifndef AVRO_LIB_GRAPHICS_APPLICATION_H_
#define AVRO_LIB_GRAPHICS_APPLICATION_H_

#include <memory>
#include <vectory>

class VertexBuffer {

};

class UniformSet {

};

class VertexArray {

	std::vector<real_t>  coordinates_; // total number of control points of the mesh : size = np
	std::vector<index_t> triangles_;   // connectivity of control points to form triangles: size = nt * (q+1)*(q+2)/2 
	std::vector<index_t> edges_;       // connectivity of control point to form edges: size = ne * (q+1)
	std::vector<real_t>  fields_;      // total number of solution controls points: size = nfield * nt * (p+1)*(p+2)/2
};

class OpenGL_VertexArray : public VertexArray {

public:
	OpenGL_VertexArray( const TopologyBase& topology ) {

	}

private:
	GLuint vertex_array_;
	GLuint coordinate_buffer_;
	GLuint triangle_buffer_;
	GLuint edge_buffer_;
	std::vector<GLuint> texture_buffers_;
};



class WebViewer_VertexArray : public VertexArray {
public:
	WebViewer_VertexArray( const TopologyBase& topology ) {

	}

	// keep track of WV id's

};

class Manager {
	virtual void write( const VertexBuffer& buffer ) = 0;
};

class OpenGL_Manager : public Manager {

	void write( const VertexBuffer& buffer );
	void draw( index_t id , const UniformSet& uniforms );

};

class WebGL_Manager : public Manager {

	void write( const VertexBuffer& buffer );
};

class Plot;

class Primitive {

public:
	Primitive( const Plot& plot );

	void create_buffers( const TopologyBase& topology ) {
		// extract triangles from topology

		// 
	}

private:
	const Plot& plot_;

	VertexArray vertex_array_;
};

class Plot {

	typedef std::shared_ptr<Primitive> Prim_ptr;

public:
	template<typename type>
	Plot( const Topology<type>& topology );

private:
	std::vector<Prim_ptr> volumes_;
	std::vector<Prim_ptr> faces_;
	std::vector<Prim_ptr> edges_;
	std::vector<Prim_ptr> nodes_;

	mat4 model_matrix_;
};

template<>
Plot::Plot( const Topology<Simplex>& topology ) {

	// check how many volumes we have
	// or check how many partitions to gather

}


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

#endif
