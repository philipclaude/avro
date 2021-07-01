#ifndef AVRO_LIB_GRAPHICS_APPLICATION_H_
#define AVRO_LIB_GRAPHICS_APPLICATION_H_

#include "avro_params.h"

#include <memory>
#include <vectory>

namespace avro
{

class UniformSet : public ParameterSet {

};

template<typename,typename> class TrianglePrimitive;

template<typename idx_t,typename flt_t>
class SolutionPrimitive {

	friend template<typename,typename> class TrianglePrimitive;

private:
	SolutionPrimitive( std::string name , coord_t order , Shader& shader ) {
		nb_basis_ = (order+1)*(order+1)/2;
		sampler_index_ = 1; // colormap is always the first sampler
	}

	void add( )

	index_t nb_basis() const { return nb_basis_; }
	coord_t order() const { return order_; }

	index_t nb() const { return solution_.size() / nb_basis_; }

	void write( idx_t texture_buffer , idx_t texture ) {

		// bind the solution values to the texture buffer
	  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , texture_buffer) );
	  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * solution_.size() , solution_.data() , GL_STATIC_DRAW) );

		// generate a texture to hold the solution buffer
		GLuint texture;
		GL_CALL( glGenTextures( 1 , &texture) );
		GL_CALL( glActiveTexture( GL_TEXTURE0 + 1 ) );
		GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture) );
		GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , texture_buffer ) );
	}

	void set_active( Shader& shader ) {
    glActiveTexture(GL_TEXTURE0 + 1 ); // colormap always comes first
		shader.setUniform( name_ , int(sampler_index_) );
	}

private:
	coord_t order_;
	index_t nb_basis_;
	idx_t buffer_;
	idx_t texture_;
	idx_t texture_index_;
	idx_t sampler_index_;
	std::string name_;

	Shader& shader_;

	std::vector<flt_t> solution_; // solution control points
};

template<typename idx_t,typename flt_t>
class TrianglePrimitive {
public:

	TrianglePrimitive( Plot& plot , coord_t order ) :
	 	order_(order),
		active_solution_("") {
		index_t nb_basis_ = (order+1)*(order+2)/2;
	}

	~TrianglePrimitive() {
		// free the buffers
		glDeleteBuffers(1,&element_buffer);
		glDeleteBuffers(1,&texture_buffer_);
		glDeleteTextures(1,&texture_);
	}

	void add_triangle( const index_t* p , index_t np );

	void write() {
		// make sure we are bound to the correct vao

		GL_CALL( glBindVertexArray( plot.vao() ) );
	}

	void render( const std::string& solution_name ) {

		if (!visible_) return;

		// use the shader
		shader_.use();

		// set the appropriate
		std::map<std::string,std::shared_ptr<SolutionPrimitive>> it;
		it = solution_.find(solution_name);
		if (it != solution_.end()) {

			// tell the shader we want to plot a solution
			shader.setUniform("with_solution",1);

			// retrieve the solution primitive
			SolutionPrimitive& solution = *it->second.get();

			// if the solution is not currently active, we need to overwrite the data in the buffer
			if (solution_name != active_solution_) solution.write(texture_buffer_,texture);

			// set the solution texture to be active in the shader
			solution.set_active(shader_);
		}

		// set any uniforms for this patch of triangles
		// for example, maybe we want geometry triangles to glow, etc.

		// bind the buffer associated with the triangle patch
		GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , element_buffer_ ) );

		// set the appropriate patch parameter
		if (order_ > 1) {
			GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );

			// draw the elements
			GL_CALL( glDrawElements( GL_PATCHES , indices_.size() ) );
		}
		else {
			GL_CALL( glDrawElements( GL_TRIANGLES , indices_.size() ) );
		}
	}

private:
	Plot& plot_;
	index_t order_;
	index_t nb_basis_;
	std::vector<idx_t> indices_;
	std::map<std::string,SolutionPrimitive> solution_;
	std::string active_solution_;

	idx_t element_buffer_;
	idx_t texture_buffer_;
	idx_t texture_;
	bool  visible_;
};

template<typename idx_t>
class EdgePrimitive {
public:
	EdgePrimitive( coord_t qorder );

	void write();
	void render();

	std::vector<idx_t> indices_;
	idx_t buffer_;
	bool  visible_;
};

template<typename idx_t,typename flt_t>
class PointPrimitive {

public:

	PointPrimitive( const Points& points );

	void write();
	void render() {
		// draw the points
	}

	idx_t buffer() const { return buffer_; }

private:
	std::vector<flt_t> coordinates_;
	idx_t buffer_;
	bool  visible_;
};

class OpenGL_VertexArray {

public:
	OpenGL_VertexArray( const TopologyBase& topology ) :
		points(topology.points()) {

		const Fields& fields = topology.fields();

		// compute all the triangle facets of the topology

		// split up the facets into triangle primitives

		// define a solution texture on each triangle patch



		// compute all the edges of the topology

		// split up the edges into edge primitives
	}

	void render() {

		// bind the vertex array
		// assign the coordinate buffer to the vertex attribute
		GL_CALL( glBindVertexArray(vertex_array_) );
		GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_.buffer()  ) );
		GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
		GL_CALL( glEnableVertexAttribArray(0) );
		GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

		// set the matrices


		// render the triangle primitives
		for (index_t k = 0; k < triangles_.size(); k++) {
			triangles_[k].render();
		}

		// loop through the edge edge_buffers
		for (index_t k = 0; k < edge_buffers_.size(); k++) {
			edges_[k].render();
		}

		// do we want to render points?
		if (points_visible_) {
			points_.render();
		}


	}

private:
	typedef GLuint  idx_t;
	typedef GLfloat flt_t;

	idx_t vertex_array_;

	PointPrimitive<idx_t,flt_t> points_;
	std::vector<TrianglePrimitive<idx_t,flt_t>> triangles_;
	std::vector<EdgePrimitive<idx_t>> edges_;
};



class WebViewer_VertexArray {
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

} // avro

#endif
