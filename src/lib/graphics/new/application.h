#ifndef AVRO_LIB_GRAPHICS_APPLICATION_H_
#define AVRO_LIB_GRAPHICS_APPLICATION_H_

class VertexBuffer {

public:
	VertexBuffer();

private:
	std::vector<real_t>  data_;
	std::vector<index_t> indices_;
};

class UniformSet {

};

class Manager {
	virtual void write( const VertexBuffer& buffer ) = 0;
};

class OpenGL_Manager : public Manager {

	void write( const VertexBuffer& buffer );
	void draw( index_t id , const UniformSet& uniforms );

	std::map<int,VertexBuffer*> buffer_;
};

class WebGL_Manager : public Manager {

	void write( const VertexBuffer& buffer );
};

class Primitive {

public:
	Primitive( const TopologyBase& topology , const Primitive* parent=nullptr );

private:
	const TopologyBase& topology_;
	const Primitive* parent_;
	mat4 model_matrix_;

	std::vector<VertexBuffer> buffers_;
	vec3 center_;
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

#endif
